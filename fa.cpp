/*
 *    Description:
 *
 *        Version:  1.0
 *        Created:  16-06-14 21:17:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Roel Kluin,
 */
#include <stdlib.h> // realloc()
#include <ctype.h> // isspace()
#include <stdio.h> // fprintf(), fseek();
#include <errno.h> // ENOMEM
#include <string.h> // memset()
#include <limits.h> //INT_MAX
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
//#include <glib.h>
#include <pcre.h>
#include "fa.h"

//required format: >ID SEQTYPE:IDTYPE LOCATION [META] (TODO: use .fai istead)
static int
parse_header(kct_t* kc, void* g, int (*gc) (void*), uint32_t b2pos)
{
    pcre *reCompiled;
    pcre_extra *pcreExtra;
    const char *str;
    unsigned start = kc->hdr_l;
    const char *re = "^(?:[^ ]+) dna:[^ :]+ [^ :]+:[^:]+:[^:]+:(\\d+):(?:\\d+):(?:\\d+)(?: .*)?$";
    int x, end = 0, subStrVec[6];
    do {
        x = gc(g);
        if (x == -1) return -2;
        _buf_grow(kc->hdr, 1ul);
        if (end == 0 && isspace(x)) end = kc->hdr_l;
        kc->hdr[kc->hdr_l++] = x;
        fputc(x, stderr);
    } while (x != '\n');
    kc->hdr[kc->hdr_l - 1] = '\0';

    reCompiled = pcre_compile(re, 0, &str, &x, NULL);
    ASSERT(reCompiled != NULL, return -3, " during Compile '%s': %s\n", re, str);

    // Optimize the regex, optional, otherwise use pcreExtra = NULL in place
    pcreExtra = pcre_study(reCompiled, 0, &str); // str is NULL unless error
    ASSERT(str == NULL, return -4, " during Study '%s': %s\n", re, str);

    x = pcre_exec(reCompiled, pcreExtra, &kc->hdr[start], kc->hdr_l - start - 1, 0, 0, subStrVec, 6);
    // x is zero on too many substrings, PCRE_ERROR_NOMATCH means no match, else error
    ASSERT(x > 0, return -5, "%d: during Extract '%s' in /%s/\n", x, &kc->hdr[start], re);

    // x - 2 is length, beide 1-based.
    pcre_get_substring(&kc->hdr[start], subStrVec, x, 1, &str);
    Bnd bd = { .b2pos = b2pos, .len = (atoi(str) - 1) | REF_CHANGE };
    //EPR("%s, %u stored", str, bd.len & NR_MASK);
    kc->bnd.push_back(bd);
    kc->hdr[end] = '\0'; // store just ID.
    kc->hdr_l = end + 1;

    pcre_free_substring(str);
    return 0;
}

/*
 * process fasta, store 2bit string and count occurances and locations of keys.
 * store per key the count and the last position.
 */
static int
fa_kc(kct_t* kc, void* g, int (*gc) (void*))
{
    uint32_t dna = 0u, rc = 0u, ndx, b2pos = 0u;
    int c;
    uint16_t t = KEY_WIDTH;
    uint8_t b = 0;
    uint8_t *s = kc->seq;
    Kcs *kcp;
    *s = '\0';

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID

    while (c >= 0) {
        b ^= b;
        switch(c) {
        case '\n': break;
        case 'C': b = 2;
        case 'G': b ^= 1;
        case 'T': b ^= 2;
        case 'A':
            dna = (dna << 2) & KEYNT_MASK; // do not set carriage bit.
            dna |= b;
            rc = ((b ^ 2)<< KEYNT_TOP) | (rc >> 2);

            *s |= b << ((b2pos & 3) << 1);
            if ((++b2pos & 3) == 0) {
                _buf_grow2(kc->seq, 1ul, s);
                *++s = '\0';
                kc->seq_l++;
            }
            if (t) {
                --t;
                break;
            }
            // ndx (key) is 2nd cNt bit excised, rc or dna
            ndx = __get_next_ndx(ndx, dna, rc);

            if (kc->kcsndx[ndx] == ~0u) {
                //EPR("%x", ndx);
                _buf_grow(kc->kcs, 1ul);
                kc->kcs[kc->kcs_l].infior = 0;
                kc->kcsndx[ndx] = kc->kcs_l++;
            }
            kcp = kc->kcs + kc->kcsndx[ndx];

            kcp->b2pos = b2pos; // replace with last position OR could keep it at first
            kcp->infior++;
            break;
        default:
            t = KEY_WIDTH; // need to rebuild key
            if (c != '>') {
                ndx = 0; // use ndx to count stretch (key is rebuilt later)
                do {
                    if (c != '\n') {
                        if (c != 'N') {
                            if (B6(b, c) || c == '>') break;
                            EPR("Ignoring strange nucleotide:%c\n", c);
                        }
                        ++ndx;
                    }
                    c = gc(g);
                } while(c != -1);
                Bnd bd = { .b2pos = b2pos, .len = ndx | N_STRETCH };
                kc->bnd.push_back(bd);
                continue; /* reprocess c */
            }
            c = parse_header(kc, g, gc, b2pos);
            if (c < 0) continue; /* reprocess c - break out of loop */
        }
        c = gc(g);
    }
    /* After this function kc->seq_l has a different function: tot no 2bits instead
     * of characters in seq, which can be derived from it: it is (kc->seq_l >> 2). */
    kc->seq_l = b2pos;
    return kc->hdr[kc->hdr_l - 1] == '\0'; // make sure last header is complete
}

/*
 * Decrement counts for keys that occur within range before a unique sequence key.
 * If another key hence has become unique, the range is extended.
 * Positions of keys with lower infior are to be considered before keys with higher
 * infior, when trying to find a position for the keys in a read.
 *
 * Returned zero if the range adjoins the downstream boundary, otherwise the start
 * position of the range.
 */
uint32_t handle_before_unique(kct_t* kc, uint32_t b2pos, uint32_t infior,
        uint32_t* fk, unsigned rest, unsigned ext)
{
    unsigned lb = rest > ext ? rest - ext : 0;
    unsigned inf_range = lb;

    while (rest != lb) {
        --rest;
        --b2pos;
        uint32_t ndx = fk[rest];
        ASSERT(ndx < kc->kcs_l, return -1u, " 0x%x", ndx);

        Kcs *z = &kc->kcs[ndx];

        // when the original boundary is reached increase the infior
        if (rest < inf_range) {
            ++infior;
            inf_range = lb;
        }
        ASSERT(z->infior != 0u, return -1u, ":Unencountered key: 0x%x", ndx);
        ASSERT((z->infior & INFERIORITY_BIT) == 0, return -1u, ":Already handled");

        if (--z->infior != 0u) { // extend if another key hence has come unique
            //EPR("%u\tndx:0x%x\tct:%u (decrementing)", b2pos, ndx, z->infior);
        } else {
            // extend if another key hence has come unique
            z->b2pos = b2pos;
            z->infior = infior;
            // alternatively increment infior for each extension; then also remove
            // z->infior = ++infior | INFERIORITY_BIT; // the rest < inf_range branch.

            lb = rest > ext ? rest - ext : 0;
        }
    }
    return rest != 0 ? b2pos : 0u;
}

/*
 * A short sequence that occurs only at one position on the genome is a `unique key'.
 * Such unique keys relate only to single positions and therefore provide a primary
 * sequence to position translation. For other keys within read- minus keylength range
 * of such hallmark sequences, it is not strictly necessary to consider this location.
 * The primary unique key should suffice to identify reads containing these keys.
 * This can lead to keys within range to become secondarily unique for positions
 * elsewhere on the genome, and so on. As long as we keep track of which key is primary
 * unique, and which secondary, or the inferiority per key, it is possible to
 * map the reads within reach of such keys.
 *
 * In this function we iterate over a genome to extend the range of unique positions
 * while keeping track of the priecedence.
 */
int extd_uniq(kct_t* kc, uint32_t* fk, unsigned gi, unsigned ext)
{
    unsigned uq = 0, fk_m, fk_l;

    uint32_t b2pos = 0u; // twobit nr
    char* hdr = kc->hdr;
    std::list<Bnd>::iterator it = kc->bnd.begin();
    int32_t addpos = 0u; // b2pos + addpos = real genomic coordinate (per chromo)

    if (gi == 0)
        fk = _buf_init(fk, 16); //init former_keys on first iteration

    while (b2pos != kc->seq_l) {
        EPR("%s:%u,%u, <start>", hdr, b2pos + addpos, it->b2pos);
        addpos += it->len & NR_MASK;
        ++it; // consider next

        uint32_t ubound = it != kc->bnd.end() ? it->b2pos : kc->seq_l;
        if (b2pos + KEY_WIDTH < ubound) { // rebuild key if we can
            EPR("%s:%u Next boundary at %u, 0x%x", hdr, b2pos + addpos, ubound, it->len);

            uint32_t i, b, dna, rc;
            fk_l = 0; // flush buffer
            __init_key(i, b, kc->seq, b2pos, dna, rc);

            // next start processing keys
            uint32_t start = ~0u;
            while (b2pos != ubound) {
                ASSERT(b2pos < ubound, return -1, "beyond ubound");
                uint32_t ndx = __next_ndx_with_b2(b, kc->seq, b2pos, dna, rc);
                //EPR("%s:%u\t0x%x => ?0x%x?", hdr, b2pos + addpos, dna, dna << 2);

                ASSERT(ndx < (1u << kc->kcsndx_m), return -1);
                ndx = kc->kcsndx[ndx];
                ASSERT(ndx < kc->kcs_l, return -1, " 0x%x", ndx);
                Kcs *z = &kc->kcs[ndx];
                /* if zero, the key was already handled in a previous iteration,
                 * but then a boundary should have been inserted to skip that
                 * region.
                 */
                ASSERT((z->infior & TYPE_MASK) == 0, return -1,
                        "key 0x%x already handled in a previous iteration", ndx);

                if (z->infior > 1) {
                    _buf_grow(fk, 1ul);
                    fk[fk_l++] = ndx;
                    ++b2pos;
                    continue;
                }
                // we encounter a unique key
                ++uq;

                z->infior = gi | INFERIORITY_BIT | // set handled, and
                    (-(dna & KEYNT_STRAND) & ORIENTATION); // central Nt state

                /* Add a new range. FIXME: if we're extending a range in a subsequent
                 * iteration the inferiority should be derived form the extended region
                 */
                unsigned infior = gi + 1u;
                start = handle_before_unique(kc, b2pos, infior, fk, fk_l, ext);
                fk_l = 0;
                ASSERT(start != -1u, return -1, "...");

                // handle after unique
                unsigned ub = min(b2pos + ext, ubound);
                unsigned inf_range = ub;
                while (++b2pos != ub) {
                    ndx = __next_ndx_with_b2(b, kc->seq, b2pos, dna, rc);

                    ASSERT(ndx < (1u << kc->kcsndx_m), return -1);
                    ASSERT(kc->kcsndx[ndx] < kc->kcs_l, return -1,
                            " kc->kcsndx[0x%x]: 0x%x", ndx, kc->kcsndx[ndx]);
                    Kcs *z = &kc->kcs[kc->kcsndx[ndx]];

                    ASSERT(z->infior != 0u, return -1, ":Already decremented");
                    ASSERT((z->infior & INFERIORITY_BIT) == 0, return -1,
                            ":Already handled %u\tdna:0x%x\tndx:0x%x\tct:%u",
                            b2pos, dna, ndx, z->infior);

                    if (b2pos >= inf_range) {
                        ++infior;
                        inf_range = ub;
                    }

                    if (--z->infior == 0u) { // extend
                        z->b2pos = b2pos;
                        ub = min(b2pos + ext, ubound);
                        //EPR("%u\tdna:0x%x\tndx:0x%x\tct:%u\t<extending til %u>",
                        //        b2pos, dna, ndx, z->infior, ub);
                        z->infior = infior | INFERIORITY_BIT;
                    }
                }
                if (b2pos != ubound || (it->len >= RESERVED)) {
                    //EPR("%u\t0x%x", start, (--it)++->len);
                    if (start != 0 || (--it)++->len >= RESERVED) { //decr safe: begin is chromo
                            //EPR("not joining either: insert (is before current it)");
                        Bnd bd = {.b2pos = start, .len = b2pos - start};
                        //EPR("%s\t%u\t%u\t(ins)", hdr, start + addpos, b2pos + addpos);
                        kc->bnd.insert(it, bd); // insert before, it stays at next
                    } else {
                        //EPR("%s\t%u\t%u\t(extd)", hdr, start + addpos, b2pos + addpos);
                        //EPR("only left joining: update b2pos and len");
                        (--it)++->len = b2pos - ub;
                        ASSERT(it != kc->bnd.begin(), return -1);
                    }
                } else {
                    EPR("right joining");
                    if (start == 0 && (--it)++->len < RESERVED) {
                        //EPR("joining both: remove one and update other");
                        start = (--it)->b2pos;
                        it = kc->bnd.erase(it); // also shifts right
                    } // else only right joining: update b2pos and len
                    ASSERT(it != kc->bnd.begin(), return -1);
                    it->len = it->b2pos + it->len - start;
                    it->b2pos = start;
                }
                //EPR("b2pos:%u\tub:%u\tubound:%u", b2pos, ub, ubound);
            }
            ASSERT(start != ~0u, return -1, "missing sequence for %s?", hdr);
        } else {
            EPR("Skipping %u\tat %u", it->b2pos, b2pos + addpos);
        }
        if (b2pos == kc->seq_l) {
            EPR("== %u uq regions in iteration ==", uq);
            if (uq == 0u)
                break;
            b2pos = addpos = 0u;
            hdr = kc->hdr;
            it = kc->bnd.begin();
            continue;
        }
        if ((it->len & TYPE_MASK) == REF_CHANGE) {
            EPR("ubound reached(2): b2pos:%u, %u, %d, %u", b2pos, kc->seq_l, addpos, ubound);
            addpos = -b2pos;
            while (*hdr++ != '\0') {} // next chromo
        }
        EPR("%s:%u,%u, <end>", hdr, b2pos + addpos, it->b2pos);
    }
    return uq;
}

// process keys and occurance in each range
// returns zero on success, negative on error
int get_all_unique(kct_t* kc, unsigned readlength)
{
    int uq, gi = 0;
    uint32_t* fk = NULL; // former_keys alloc on first iteration in extd_uniq();
    do {
        uq = extd_uniq(kc, fk, gi, readlength - KEY_WIDTH);
        EPR("%i unique regions extended in genome iteration %u", uq, ++gi);
    } while (uq > 0);
    free(fk);
    return uq;
}

int
fa_print(seqb2_t *fa)
{
    return 1;
}

int
fn_convert(struct gzfh_t* fhout, const char* search, const char* replace)
{
    char* f = strstr(fhout->name, search);
    if (f == NULL) return 0;
    strncpy(f, replace, strlen(replace) + 1);
    return 1;
}

int fa_index(struct seqb2_t* seq)
{
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    struct gzfh_t* fhin = seq->fh + 2; // init with ref
    uint64_t blocksize = (uint64_t)seq->blocksize << 20;
    int res, ret = -1;
    void* g;
    int (*gc) (void*);
    int is_gzfile = fhin->io != NULL;
    char file[256];

    const char* ndxact[4] = {".fa", "_kcpos.wig.gz"};

    if (fhout->name != NULL && ((res = strlen(fhout->name)) > 255)) {
        strncpy(file, fhout->name, res);
    } else {
        res = strlen(fhin->name);
        ASSERT(strncmp(fhin->name, "stdin", res) != 0, return -1);
        strcpy(file, fhin->name);
    }
    fhout->name = &file[0];

    kct_t kc;
    kc.H = kh_init(UQCT);
    kc.seq = _buf_init(kc.seq, 16);
    kc.hdr = _buf_init(kc.hdr, 8);
    kc.kcs = _buf_init(kc.kcs, 16);
    kc.kcsndx = _buf_init_arr(kc.kcsndx, KEYNT_BUFSZ_SHFT);
    // the first entry (0) is for a reference ID always.
    memset(kc.kcsndx, ~0u, KEYNT_BUFSZ * sizeof(kc.kcsndx[0]));

    if (is_gzfile) {
        g = fhin->io;
        gc = (int (*)(void*))&gzgetc;
    } else {
        g = fhin->fp;
        gc = (int (*)(void*))&fgetc;
    }
    /* TODO: load dbSNP and known sites, and mark them. */
    res = fa_kc(&kc, g, gc);
    ASSERT(res >= 0, goto out);
    res = get_all_unique(&kc, seq->readlength);
    ASSERT(res >= 0, goto out);

    ret = res != 0;

out:
    _buf_free(kc.kcs);
    _buf_free(kc.kcsndx);
    _buf_free(kc.hdr);
    _buf_free(kc.seq);
    kh_destroy(UQCT, kc.H);
    return ret;

}


