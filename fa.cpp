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
    uint32_t dna = 0u, rc = 0u;
    uint32_t ndx, b2pos = 0u, b = 0ul;
    int c;
    uint16_t t = KEY_WIDTH;
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

            *s |= b << ((b2pos & 3) << 1);
            dna = _seq_next(b, dna, rc);
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
            ndx = _get_ndx(ndx, dna, rc);

            if (kc->kcsndx[ndx] == UNINITIALIZED) {
                //EPR("%x", ndx);
                _buf_grow(kc->kcs, 1ul);
                kc->kcs[kc->kcs_l].infior = 0;
                kc->kcs[kc->kcs_l].ct = 0;
                kc->kcsndx[ndx] = kc->kcs_l++;
            }
            kcp = kc->kcs + kc->kcsndx[ndx];

            kcp->b2pos = b2pos; // replace with last position OR could keep it at first
            kcp->ct++;
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
    Bnd bd = { .b2pos = b2pos, .len = -1u };
    //EPR("%s, %u stored", str, bd.len & NR_MASK);
    kc->bnd.push_back(bd);
    return kc->hdr[kc->hdr_l - 1] == '\0'; // make sure last header is complete
}

// XXX: problem if w->ct != 1: w->b2pos may refer to later
// b2pos, with a different 'Next b2', solutions:
// look to next seq, look further previous, untill we can confirm it?
// walk all sequences again?
uint32_t last_b2pos_for_ndx(kct_t* kc, uint32_t odna, uint32_t orc, uint32_t oct,
        uint32_t b)
{
    uint32_t dna, rc, i, Nt, b2pos = 0, ct = oct;
    for (Nt = 0; Nt != 4; ++Nt) {
        dna = odna;
        rc = orc;
        dna = _seq_prev(Nt, dna, rc);
        i = _get_ndx(i, dna, rc);
        if (kc->kcsndx[i] == UNINITIALIZED) continue;

        Kcs* w = &kc->kcs[kc->kcsndx[i]];
        if (w->ct != 1u) break;
        i = w->b2pos + 1;
        if (read_b2(kc->seq, i) == b) {
            if (i > b2pos) b2pos = i;
            if (--ct == 0) return b2pos;
        }
    }
    ct = oct;
    for (Nt = 0; Nt != 4; ++Nt) {
        dna = odna;
        rc = orc;
        dna = _seq_next(Nt, dna, rc);
        i = _get_ndx(i, dna, rc);
        if (kc->kcsndx[i] == UNINITIALIZED) continue;

        Kcs* w = &kc->kcs[kc->kcsndx[i]];
        if (w->ct != 1u) break;
        i = w->b2pos - 1;
        if (read_b2(kc->seq, i) == b) {
            if (i > b2pos) b2pos = i;
            if (--ct == 0) return b2pos;
        }
    }
    return UNINITIALIZED;
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
 * while keeping track of the precedence. it returns the number of keys that have now
 * become unique.
 */


int extd_uniq(kct_t* kc, uint32_t* fk, unsigned ext)
{
    unsigned i, ct = 0, dbg = 0;
    Kcs *x = NULL, *y;
    const Kcs *kcs_end = &kc->kcs[kc->kcs_l];
    std::queue<uint32_t> updatepos;

    for (y = kc->kcs; y != kcs_end; ++y) {
        if (y->ct != 1u)
            continue;
        if (x == NULL) { // first
            x = y;
            continue;
        }
        // for unique positions, the b2pos must be correct. This is true in the
        // first iteration
        if ((y->b2pos - x->b2pos) < ext) {
            uint32_t dna = 0u, rc = 0u;
            uint32_t b2pos = x->b2pos - KEY_WIDTH, b = 0u;
            _init_key(i, b, kc->seq, b2pos, dna, rc);
            uint32_t end = y->b2pos - 1;

            while (b2pos++ != end) {
                uint32_t ndx = _read_next_ndx(b, kc->seq, b2pos, dna, rc);
                Kcs* z = &kc->kcs[kc->kcsndx[ndx]];
                if (--(z->ct) == 1u)
                    ++ct;
                if (b2pos != z->b2pos)
                    continue;
                // the last b2pos (registered for ndx) was removed. Search next last
                // b2pos for ndx
                i = last_b2pos_for_ndx(kc, dna, rc, z->ct, b);
                if (i != UNINITIALIZED) {
                    z->b2pos = i;
                    continue;
                }
            }
        } else {
            // insert range
            Rng rg = { .s = x->b2pos, .e = y->b2pos };
            kc->rng.push_back(rg);
            x = NULL;
        }
    }
    // push last range

    // next iterate back
    i = (x->b2pos - KEY_WIDTH - 1) ^ NEEDS_UPDATE;
    while(ct) {
    }
    
    // second iteration: walk all sequences, except already processed ranges
    // and restore positions.
}

// process keys and occurance in each range
// returns zero on success, negative on error
int get_all_unique(kct_t* kc, unsigned readlength)
{
    int uq, gi = 0;
    uint32_t* fk = NULL; // former_keys alloc on first iteration in extd_uniq();
    do {
        uq = extd_uniq(kc, fk, readlength - KEY_WIDTH);
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
    memset(kc.kcsndx, UNINITIALIZED, KEYNT_BUFSZ * sizeof(kc.kcsndx[0]));

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


