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

enum ensembl_parts {ID,  SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META};

enum bnd_type {NEW_REF_ID, N_STRETCH, UQ_REGION, END_REF, DBSNP, KNOWN_SITE};

//ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
// fai does not handle chromosomes with offset.
Hdr *new_header(kct_t* kc, void* g, int (*gc) (void*), Bnd *bd)
{
    int c;
    Hdr* h = new Hdr;
    h->hdr_type = ENSEMBL_HDR;

    h->id = _buf_init_err(h->id, 2, return NULL);
    unsigned p = ID;
    h->part[p] = h->id;
    while ((c = gc(g)) != -1) {
        fputc(c, stderr);
        switch (c) {
            case '\n': break;
            case ' ':
                if (p == IDTYPE) {
                    char* q = h->part[p];
                    if (strncmp(q, "chromosome", strlen(q)) != 0 &&
                            strncmp(q, "supercontig", strlen(q)) != 0 &&
                            strncmp(q, "nonchromosomal", strlen(q)) != 0)
                        h->hdr_type = UNKNOWN_HDR;
                } else if (p != ID) {
                    h->hdr_type = UNKNOWN_HDR;
                }
                h->id[h->id_l++] = '\0';
                if (p < ARRAY_SIZE(h->part))
                    h->part[++p] = h->id;
                continue;
            case ':': 
                if (p == ID || p == IDTYPE)
                    h->hdr_type = UNKNOWN_HDR;
                else if (p == SEQTYPE) {
                    char* q = h->part[p];
                    if (strncmp(q, "dna", strlen(q)) != 0) {
                        if (strncmp(q, "dna_rm", strlen(q)) != 0) {
                            fprintf(stderr, "\nWARNING: reference isrepeatmasked.\n");
                        } else if (strncmp(q, "cdna", strlen(q)) != 0) {
                            fprintf(stderr, "\nWARNING: aligning against cDNA.\n");
                        } else if (strncmp(q, "pep", strlen(q)) != 0) {
                            fprintf(stderr, "\nERROR: peptide alignment?\n");
                            delete h;
                            return NULL;
                        } else if (strncmp(q, "rna", strlen(q)) != 0) {
                            fprintf(stderr, "\nWARNING: aligning against rna.\n");
                        } else { // could still be ok.
                            h->hdr_type = UNKNOWN_HDR;
                        }
                    }
                }
                _buf_grow_err(h->id, 1ul, return NULL);
                h->id[h->id_l++] = '\0';
                if (p < ARRAY_SIZE(h->part))
                    h->part[++p] = h->id;
                continue;
            default:
                if (isspace(c)) {
                    c = '\0';
                    h->hdr_type = UNKNOWN_HDR;
                }
                _buf_grow_err(h->id, 1ul, return NULL);
                h->id[h->id_l++] = c;
                continue;
        }
        break;
    }
    if (p < NR) 
        h->hdr_type = UNKNOWN_HDR;

    _buf_grow_err(kc->bd, 1ul, return NULL);
    bd = &kc->bd[kc->bd_l++];
    ASSIGN_BD(bd, NEW_REF_ID, ~0u, 0, ~0u, 0, 0);

    if (h->hdr_type == ENSEMBL_HDR) {
        h->p_l = p;
        p = atoi(h->part[START]);
        if (p != 1)
            bd->s = p;
        bd->l = atoi(h->part[END]);
    } else {
        h->p_l = 1;
        fprintf(stderr, "\nWARNING: non-ensembl reference.\n");
    }
    kc->hdr.insert(Map::value_type(h->id, h));
    return h;
}

void free_kc(kct_t* kc)
{
    for (Kct* k = kc->kct; k != &kc->kct[kc->kct_l]; ++k)
        if (k->seq.m > 3)
            free(k->p.b2);

    for (Hdr* h = kc->h; h != &kc->h[kc->h_l]; ++h)
        free(h->id);

    _buf_free(kc->h);
    _buf_free(kc->bd);
    _buf_free(kc->kct);
    _buf_free(kc->kcsndx);
}

/*
 * proh.p_lcess fasta, store 2bit string and count occurances and locations of keys.
 * store per key the count and the last position.
 */
static int
fa_kc(kct_t* kc, void* g, int (*gc) (void*))
{
    uint64_t dna = 0, rc = 0, ndx, b = 0;
    uint32_t b2pos = 0u, addpos = 0u, endpos = 0u;
    int c;
    uint16_t t = KEY_WIDTH;
    Kct *ct = NULL;
    uint8_t *q = NULL;
    Hdr* h = NULL;
    Bnd* bd = NULL;

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID

    // TODO: realpos based dbsnp or known site boundary insertion.
    while (c >= 0) {
        DEBUG_ASSIGN_ENDDNA(bd, dna);
        if (c != '>') { // N-stretch
            t = 0; // use  to count stretch (key is rebuilt later)
            do {
                if (c != '\n') {
                    if (c != 'N') {
                        if (B6(b, c) || c == '>') break;
                        EPR("Ignoring strange nucleotide:%c\n", c);
                    }
                    ++t;
                }
                c = gc(g);
            } while(c != -1);
            _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
            bd = &kc->bd[kc->bd_l++];
            ASSIGN_BD(bd, N_STRETCH, ~0u, b2pos + addpos, t, 0, dna);
            addpos += t;
        } else { // header
            if (h != NULL) {
                _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
                bd = &kc->bd[kc->bd_l++];
                ASSIGN_BD(bd, END_REF, ~0u, b2pos + addpos, 0, 0, dna);
                h->bnd.push_back(bd);
            }
            b2pos = 0u;
            Hdr* h = new_header(kc, g, gc, bd);
            if (h == NULL) return -ENOMEM;
            addpos = bd->s;
            endpos = bd->l;
            c = gc(g);
        }
        for (t = b2pos + KEY_WIDTH; b2pos != t; ++b2pos) {
            b ^= b;
            switch(c) {
                case 'C': case 'c': b = 2;
                case 'G': case 'g': b ^= 1;
                case 'U': case 'u': case 'T': case 't': b ^= 2;
                case 'A': case 'a': dna = _seq_next(b, dna, rc);
                case '\n':
                    c = gc(g);
                    continue;
            }
            break;
        }
        if_ever (b2pos != t) // c == -1, 'N' or '>'
            continue;
        // ndx (key) is 2nd cNt bit excised, rc or dna
        bd->dna = dna;
        h->bnd.push_back(bd);
        ndx = _get_ndx(ndx, b, dna, rc);
        do {
            if (kc->kcsndx[ndx] != UNINITIALIZED) {
                ct = kc->kct + kc->kcsndx[ndx];
                if (ct->seq.m == 3) { // Kct has seq format;
                    q = &ct->seq.b2[ct->seq.l >> 2];
                    if ((ct->seq.l++ & 3) != 0) {
                        *q <<= 2;
                    } else {
                        b = ARRAY_SIZE(ct->seq.b2);
                        if (ct->seq.l == (b << 2)) { // convert to p (pointer) format
                            q = (uint8_t *)malloc(16);
                            memcpy(q, ct->seq.b2, b); 
                            ASSERT(ct->seq.m == ct->p.m, return -EINVAL,
                                    "struct not aligned as expected\n");
                            ++ct->p.m;
                            ASSERT(ct->seq.l == (ct->p.l & 0xff), return -EINVAL,
                                    "redundant endian test failed\n");
                            ct->p.l = b << 2;
                            ct->p.b2 = q;
                            q += b;
                        }
                        *q = '\0';
                    }
                } else { // Kct in p format;
                    q = ct->p.b2 + (ct->p.l >> 2);
                    if ((ct->p.l++ & 3) != 0) {
                        *q <<= 2;
                    } else {
                        if (ct->p.l == (1u << ct->p.m)) {
                            q = (uint8_t *)realloc(ct->p.b2, 1 << ++ct->p.m);
                            if_ever (q == NULL) return -ENOMEM;
                            ct->p.b2 = q;
                        }
                        *q = '\0';
                    }
                }
            } else { // new key
                //EPR("%x", ndx);
                _buf_grow(kc->kct, 1ul);
                kc->kcsndx[ndx] = kc->kct_l;
                ct = kc->kct + kc->kct_l++;
                ct->seq = {.m = 3, .l = 1, 0}; // init remaining to zero
            }
            // next iteration ...
            while (isspace(c = gc(g)));
            b ^= b;
            switch(c) {
                case 'C': case 'c': b = 2;
                case 'G': case 'g': b ^= 1;
                case 'T': case 'U': case 't': case 'u': b ^= 2;
                    *q |= b; // append next 2bit to last ndx
                case 'A': case 'a': dna = _seq_next(b, dna, rc);
            }
        } while (c >= 0);
        ASSERT(b2pos + addpos < endpos, return -ERANGE);
    }
    if (h != NULL) {
        DEBUG_ASSIGN_ENDDNA(bd, dna);
        _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
        bd = &kc->bd[kc->bd_l++];
        ASSIGN_BD(bd, END_REF, ~0u, b2pos + addpos, 0, 0, dna);
        h->bnd.push_back(bd);
    }
    //EPR("%s, %u stored", str, bd.len & NR_MASK);
    return 0;
}

/* Reverse Complement */
static inline uint64_t revcmp(uint64_t dna)
{
    uint64_t m = 0x3333333333333333UL;
    dna = ((dna & m) << 2) | ((dna & (~m)) >> 2);
    m = 0x0f0f0f0f0f0f0f0fUL;
    dna = ((dna & m) << 4) | ((dna & (~m)) >> 4);
    asm ("bswap %0" : "=r" (dna) : "0" (dna));
    dna ^= 0xaaaaaaaaaaaaaaaa;
    return dna >> (64 - (KEY_WIDTH << 1));
}


int extd_uniq(kct_t* kc, uint32_t* fk, unsigned ext)
{
    uint64_t t = kc->kct_l * sizeof(Walker), dna = 0ul;
    Walker* w = (Walker*)malloc(t);
    memset(w, 0ul, t);
    uint32_t uqct, iter = 0;

    do {
        uqct = 0;
        for (Hdr* h = kc->h; h != &kc->h[kc->h_l]; ++h) {
            std::list<Bnd*>::iterator bd = h->bnd.begin();


            uint32_t pos = (*bd)->s, endpos = (*bd)->l;
            do {
                ASSERT(pos < endpos, return -EINVAL);
                dna = (*bd)->dna;
                uint64_t rc = revcmp(dna);
                uint32_t remain = ext;
                uint32_t infior = 1;
                Kct *x;

                for (;pos != (*bd)->s; ++pos) {
                    uint64_t ndx = _get_ndx(ndx, t, dna, rc);

                    Kct *y = kc->kct + kc->kcsndx[ndx];
                    if (remain) {
                        if (--remain) {
                            if (infior > w[ndx].infior)
                                w[ndx].infior = infior;
                        } else {
                            infior = 0;
                        }
                    }
                    t = w[ndx].count;
                    if (y->seq.m > 3) {
                        ASSERT(t < y->p.l, return -EINVAL);
                        t = y->p.b2[t >> 2] >> (t << 1);
                        w[ndx].count++;
                        continue;
                    }
                    ASSERT(t < y->seq.l, return -EINVAL);
                    t = y->seq.b2[t >> 2] >> (t << 1);
                    if (y->seq.l > 1) {
                        w[ndx].count++;
                        continue;
                    }
                    // found unique
                    if (remain == 0) {
                        // found first unique
                        x = y;
                        remain = ext;
                        infior = w[ndx].infior + 1;
                        continue;
                    }
                    ++uqct;
                    // found second unique
                }
#ifdef TEST_BND
                ASSERT((*bd)->end_dna = dna, return -EFAULT);
#endif
                pos += (*bd)->l;
            } while (++bd != h->bnd.end());
            EPR("extended %u uinique ranges in iteration %u", uqct, ++iter);
        }
    } while (uqct != 0);
    free(w);
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
    kc.kct = _buf_init(kc.kct, 16);
    kc.h = _buf_init(kc.h, 0);
    kc.bd = _buf_init(kc.bd, 0);
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
    free_kc(&kc);
    return ret;

}


