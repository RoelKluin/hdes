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

static const uint8_t u8revcmp[256] = { // TODO: check
    0xaa, 0xea, 0x2a, 0x6a, 0xba, 0xfa, 0x3a, 0x7a, 0x8a, 0xca, 0x0a, 0x4a, 0x9a, 0xda, 0x1a, 0x5a,
    0xae, 0xee, 0x2e, 0x6e, 0xbe, 0xfe, 0x3e, 0x7e, 0x8e, 0xce, 0x0e, 0x4e, 0x9e, 0xde, 0x1e, 0x5e,
    0xa2, 0xe2, 0x22, 0x62, 0xb2, 0xf2, 0x32, 0x72, 0x82, 0xc2, 0x02, 0x42, 0x92, 0xd2, 0x12, 0x52,
    0xa6, 0xe6, 0x26, 0x66, 0xb6, 0xf6, 0x36, 0x76, 0x86, 0xc6, 0x06, 0x46, 0x96, 0xd6, 0x16, 0x56,
    0xab, 0xeb, 0x2b, 0x6b, 0xbb, 0xfb, 0x3b, 0x7b, 0x8b, 0xcb, 0x0b, 0x4b, 0x9b, 0xdb, 0x1b, 0x5b,
    0xaf, 0xef, 0x2f, 0x6f, 0xbf, 0xff, 0x3f, 0x7f, 0x8f, 0xcf, 0x0f, 0x4f, 0x9f, 0xdf, 0x1f, 0x5f,
    0xa3, 0xe3, 0x23, 0x63, 0xb3, 0xf3, 0x33, 0x73, 0x83, 0xc3, 0x03, 0x43, 0x93, 0xd3, 0x13, 0x53,
    0xa7, 0xe7, 0x27, 0x67, 0xb7, 0xf7, 0x37, 0x77, 0x87, 0xc7, 0x07, 0x47, 0x97, 0xd7, 0x17, 0x57,
    0xa8, 0xe8, 0x28, 0x68, 0xb8, 0xf8, 0x38, 0x78, 0x88, 0xc8, 0x08, 0x48, 0x98, 0xd8, 0x18, 0x58,
    0xac, 0xec, 0x2c, 0x6c, 0xbc, 0xfc, 0x3c, 0x7c, 0x8c, 0xcc, 0x0c, 0x4c, 0x9c, 0xdc, 0x1c, 0x5c,
    0xa0, 0xe0, 0x20, 0x60, 0xb0, 0xf0, 0x30, 0x70, 0x80, 0xc0, 0x00, 0x40, 0x90, 0xd0, 0x10, 0x50,
    0xa4, 0xe4, 0x24, 0x64, 0xb4, 0xf4, 0x34, 0x74, 0x84, 0xc4, 0x04, 0x44, 0x94, 0xd4, 0x14, 0x54,
    0xa9, 0xe9, 0x29, 0x69, 0xb9, 0xf9, 0x39, 0x79, 0x89, 0xc9, 0x09, 0x49, 0x99, 0xd9, 0x19, 0x59,
    0xad, 0xed, 0x2d, 0x6d, 0xbd, 0xfd, 0x3d, 0x7d, 0x8d, 0xcd, 0x0d, 0x4d, 0x9d, 0xdd, 0x1d, 0x5d,
    0xa1, 0xe1, 0x21, 0x61, 0xb1, 0xf1, 0x31, 0x71, 0x81, 0xc1, 0x01, 0x41, 0x91, 0xd1, 0x11, 0x51,
    0xa5, 0xe5, 0x25, 0x65, 0xb5, 0xf5, 0x35, 0x75, 0x85, 0xc5, 0x05, 0x45, 0x95, 0xd5, 0x15, 0x55
};

//>ID SEQTYPE:IDTYPE LOCATION [META]
static int
parse_header(kct_t* kc, void* g, int (*gc) (void*), uint32_t b2pos)
{
    pcre *reCompiled;
    pcre_extra *pcreExtra;
    const char *str;
    unsigned start = kc->hdr_l;
    const char *aStrRegex = "^(?:[^ ]+) dna:[^ :]+ [^ :]+:[^:]+:[^:]+:(\\d+):(?:\\d+):(?:\\d+)(?: .*)?$";
    int x, end = 0, subStrVec[6];
    do {
        x = gc(g);
        if (x == -1) return -2;
        _buf_grow(kc->hdr, 1ul);
        kc->hdr[kc->hdr_l++] = x;
        if (end == 0 && isspace(x)) end = kc->hdr_l;
        fputc(x, stderr);
    } while (x != '\n');
    kc->hdr[kc->hdr_l - 1] = '\0';

    reCompiled = pcre_compile(aStrRegex, 0, &str, &x, NULL);
    if(reCompiled == NULL) {
        printf("ERROR: Could not compile '%s': %s\n", aStrRegex, str);
        return -3;
    }
    // Optimize the regex
    pcreExtra = pcre_study(reCompiled, 0, &str);
    if(str != NULL) {
        printf("ERROR: Could not study '%s': %s\n", aStrRegex, str);
        return -4;
    }
    x = pcre_exec(reCompiled, pcreExtra, &kc->hdr[start], kc->hdr_l - start - 1, 0, 0, subStrVec, 6);
    if(x < 0) {
        if(x != PCRE_ERROR_NOMATCH) return -5; // something bad happened
        printf("Cannot extract start from fasta header\n");
        fprintf(stderr, "matching '%s' to /%s/\n", &kc->hdr[start], aStrRegex);
        return -6;
    }
    if(x == 0) {
        printf("Too many substrings were found to fit in subStrVec!\n");
        return -7;
    }
    // x - 2 is length, beide 1-based.
    pcre_get_substring(&kc->hdr[start], subStrVec, x, 1, &str);
    Bnd bd = { .b2pos = b2pos, .len = (atoi(str) - 1) | REF_CHANGE };
    kc->bnd.push_back(bd);
    kc->hdr[end] = '\0'; // store just ID.
    kc->hdr_l = end + 1;

    pcre_free_substring(str);
    return 0;
}

/*
 * process fasta, store 2bit string and count occurances and first locations of keys.
 * store per key the count and the last position.
 */
static int
fa_kc(kct_t* kc, void* g, int (*gc) (void*))
{
    uint32_t dna = 0u, rev = 0u, ndx, b2pos = 0u;
    int c;
    uint16_t t = KEY_WIDTH;
    uint8_t b = 0;
    uint8_t *s = kc->seq;
    Kcs *kcp;
    *s = '\0';

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID

    while (c >= 0) {
        switch(c) {
        case '\n': break;
        case 'C': b = 2;
        case 'G': b ^= 1;
        case 'T': b ^= 2;
        case 'A':
            dna = (dna << 2) & KEYNT_MASK; // do not set carriage bit.
            dna |= b;
            rev = (b << KEYNT_TOP) | (rev >> 2);

            *s |= b;
            if ((++b2pos & 3) == 0) {
                _buf_grow2(kc->seq, 1ul, s);
                *++s = '\0';
                kc->seq_l++;
            }
            *s <<= 2;
            if (t) {
                --t;
                break;
            }
            // ndx (key) is 2nd cNt bit excised, rev or dna
            ndx = (dna & KEYNT_STRAND) ? dna : (rev ^ 0xaaaaaaaa);
            ndx = (((ndx >> 1) & KEYNT_TRUNC_UPPER) | (ndx & HALF_KEYNT_MASK));

            if (kc->kcsndx[ndx] == ~0u) {
                //fprintf(stderr, "%x\n", ndx);
                _buf_grow(kc->kcs, 1ul);
                kc->kcs[kc->kcs_l].ct = 0;
                kc->kcsndx[ndx] = kc->kcs_l++;
            }
            kcp = kc->kcs + kc->kcsndx[ndx];

            kcp->b2pos = b2pos; // replace with last position OR could keep it at first
            kcp->ct++;
            break;
        default:
            t = KEY_WIDTH; // need to rebuild key
            if (c != '>') {
                ndx = 0; // use ndx to count stretch (key is rebuilt later)
                do {
                    if (c != '\n') {
                        if (c != 'N') {
                            if (B6(b, c) || c == '>') break;
                            fprintf(stderr, "Ignoring strange nucleotide:%c\n", c);
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
            if (c < 0) continue; // causes break out of loop
        }
        c = gc(g);
        b ^= b;
    }
    kc->kcs_l = b2pos;
    /* TODO: fix last (what?) */
    return kc->hdr[kc->hdr_l - 1] == '\0'; // make sure last header is complete
}

#define __append_next_b2(b, s, b2pos, dna, rc) ({\
    b = b2pos;\
    b = (s[b>>2] >> ((b & 3) << 1)) & 3;\
    rc = ((b ^ 2) << KEYNT_TOP) | (rc >> 2);\
    ((dna << 2) & KEYNT_MASK) | b;\
})

#define __get_next_ndx(ndx, dna, rc) ({\
        ndx = (dna & KEYNT_STRAND) ? dna : rc;\
        ((ndx >> 1) & KEYNT_TRUNC_UPPER) | (ndx & HALF_KEYNT_MASK);\
})

#define __next_ndx_with_b2(b, s, b2pos, dna, rc) ({\
    dna = __append_next_b2(b, s, b2pos, dna, rc);\
    __get_next_ndx(b, dna, rc);\
})

#define __init_key(i, b, s, b2pos, dna, rc) ({\
    dna = rc = 0; \
    for (i = b2pos + KEY_WIDTH; b2pos != i; ++b2pos) {\
        dna = __append_next_b2(b, s, b2pos, dna, rc);\
    }\
})

/*
 * Decrement counts for keys that occur within range before a unique sequence key.
 * If another key hence has become unique, the range is extended.
 * Positions of keys with lower inferiority are to be considered before keys with higher
 * inferiority, when trying to find a position for the keys in a read.
 *
 * Returned zero if the range adjoins the downstream boundary, otherwise the start
 * position of the range.
 */
uint32_t handle_before_unique(kct_t* kc, uint32_t b2pos, uint32_t inferiority,
        uint32_t* fk, unsigned rest, unsigned ext)
{
    unsigned lb = rest > ext ? rest - ext : 0;
    unsigned inf_range = lb;

    while (--rest > lb) {
        --b2pos;
        Kcs *z = &kc->kcs[kc->kcsndx[fk[rest]]];

        // when the original boundary is reached increase the inferiority
        if (rest < inf_range) {
            ++inferiority;
            inf_range = lb;
        }

        // all should not yet be handled as unique
        assert((z->ct & INFERIORITY_BIT) == 0);

        if (--z->ct == 0) { // extend if another key hence has come unique
            z->b2pos = b2pos;
            z->ct = inferiority;
            // alternatively increment inferiority for each extension; then also remove
            // z->ct = ++inferiority | INFERIORITY_BIT; // the rest < inf_range branch.

            lb = rest > ext ? rest - ext : 0;
        }
    }
    return rest != 0 ? b2pos : 0u;
}

// decrement within range after unique
/*
 * Decrement counts for keys that occur within range after a unique sequence key.
 * If another key hence has become unique, the range is extended.
 * Positions of keys with lower inferiority are to be considered before keys with higher
 * inferiority, when trying to find a position for the keys in a read.
 *
 * Returned kc->seq_l if the range adjoins the upstream boundary, otherwise the end
 * position of the range.
 */
uint32_t handle_after_unique(kct_t* kc, uint32_t ubound, unsigned ext,
        uint32_t b2pos, uint32_t inferiority, uint32_t dna, uint32_t rc)
{
    unsigned ub = min(b2pos + ext, ubound);
    unsigned inf_range = ub;

    uint8_t* s = &kc->seq[b2pos >> 2];

    unsigned b;
    while (++b2pos < ub) {
        uint32_t ndx = __next_ndx_with_b2(b, s, b2pos, dna, rc);
        Kcs *z = &kc->kcs[ndx];
        assert(z->ct != 0 && (z->ct & INFERIORITY_BIT) == 0);
        if (b2pos >= inf_range) {
            ++inferiority;
            inf_range = ub;
        }

        if (--z->ct == 1) { // extend
            z->b2pos = b2pos;
            z->ct = inferiority;
            ub = min(ub + ext, ubound);
        }
    }

    return b2pos != ub ? b2pos : ubound;
}

/*
 * A short sequence that occurs only at one position on the genome is a `unique key'.
 * Such unique keys relate only to single positions and therefore provide a primary
 * sequence to position translation. For other keys within readlength-range of such
 * hallmark sequences, their locations therefore no longer have to be considered.
 * The hallmark should already be enough to identify reads containing this key.
 * This can lead to keys that are secondarily unique, and so on.
 *
 * In this function we iterate over a genome to extend the range of unique positions
 * while keeping track of the priecedence.
 */
int extd_uniq(kct_t* kc, uint32_t* fk, unsigned gi, unsigned ext)
{
    unsigned uq = 0, fk_m, fk_l;

    uint32_t b2pos = 0u; // twobit nr
    uint8_t* s = kc->seq;
    char* hdr = kc->hdr; // TODO
    std::list<Bnd>::iterator it = kc->bnd.begin();
    uint32_t addpos = it->len & NR_MASK; // TODO: b2pos + addpos = real genomic coordinate (per chromo)

    if (gi == 0)
        fk = _buf_init(fk, 16); //init former_keys on first iteration

    while (++it != kc->bnd.end()) { // look ahead to next mark

        if (b2pos + KEY_WIDTH < it->b2pos) { // rebuild key if we can


            uint32_t i, b, dna, rc;
            fk_l = 0; // flush buffer
            __init_key(i, b, s, b2pos, dna, rc);
            addpos += KEY_WIDTH; // construction moved us.

            do {
                uint32_t ndx = __next_ndx_with_b2(b, s, b2pos, dna, rc);
                Kcs *z = &kc->kcs[kc->kcsndx[ndx]];

                if (z->ct <= 1) {
                    // we encounter a unique key
                    ++uq;

                    /* if zero, the key was already handled in a previous iteration,
                     * but then a boundary should have been inserted to skip that
                     * region.
                     */
                    assert ((z->ct & INFERIORITY_BIT) == 0);
                    z->ct = gi | INFERIORITY_BIT; // set handled.

                    // add a new range.
                    uint32_t start = handle_before_unique(kc, b2pos, z->ct + 1, fk, fk_l, ext);

                    uint32_t ubound = it != kc->bnd.end() ? it->b2pos : kc->seq_l;
                    uint32_t end = handle_after_unique(kc, ubound, ext, b2pos, z->ct + 1, dna, rc);

                    if (end != ubound) {
                        if (start != 0 || (--it)++->len >= RESERVED) { //begin is chromo
                            //not joining either: insert (is before current it)
                            it = kc->bnd.insert(it, {.b2pos = start, .len = end - start});
                        } else {
                            //only left joining: update b2pos and len
                            --it->len = end - it->b2pos;
                        }
                    } else {
                        if (start == 0 && (--it)++->len < RESERVED) {
                            // joining both: remove one and update other
                            start = --it->b2pos;
                            it = kc->bnd.erase(it); // also shifts right
                        } // else only right joining: update b2pos and len
                        it->len = it->b2pos + it->len - start;
                        it->b2pos = start;
                    }
                    fprintf(stderr, "Unique at %u\t0x%x\n", b2pos, ndx);
                    break;
                }

                _buf_grow(fk, 1ul);
                fk[fk_l++] = ndx;
            } while (++b2pos != it->b2pos);
        } else {
            // skip
            fprintf(stderr, "Skipping %u\tat %u\n", it->b2pos - b2pos, b2pos);
            b2pos = it->b2pos;
            // TODO: also ensure none of the key buffer is used
        }
        if ((it->len & ~NR_MASK) == REF_CHANGE) {
            addpos = 0u;
            while (*hdr++ != '\0') {} // next chromo
        }
        addpos += it->len & NR_MASK;
    }
    // TODO: handle last keys
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
        fprintf(stderr, "%i unique positions added in genome iteration %u\n", uq, ++gi);
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
        if (strncmp(fhin->name, "stdin", res) == 0) {
            fputs("Cannot index from stdin (seek)\n", stderr);
            return -1;
        }
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
    if (res < 0) { fprintf(stderr, "== index failed:%d", res);  goto out; }
    res = get_all_unique(&kc, seq->readlength);
    if (res < 0) { fprintf(stderr, "== counting failed:%d", res);  goto out; }

    ret = res != 0;

out:
    _buf_free(kc.kcs);
    _buf_free(kc.kcsndx);
    _buf_free(kc.hdr);
    _buf_free(kc.seq);
    kh_destroy(UQCT, kc.H);
    return ret;

}


