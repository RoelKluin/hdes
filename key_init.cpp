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
#include "fa.h"
#include <ctype.h> // isspace()

#define ENS_HDR_PARTCT 10
//enum ensembl_parts {ID, SEQTYPE, IDTYPE, IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR};
//ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
// fai does not handle chromosomes with offset.
static Hdr*
new_header(kct_t* kc, void* g, int (*gc) (void*))
{
    int c;
    Hdr* h = new Hdr;
    h->part = (uint32_t*)malloc(ENS_HDR_PARTCT * sizeof(uint32_t));
    ASSERT(h->part != NULL, return NULL);

    kc->h.push_back(h);

    uint32_t p = ID;
    h->part[p] = kc->id_l;
    while ((c = gc(g)) != -1) {
        //fputc(c, stderr);
        switch (c) {
            case '\n': _buf_grow_add_err(kc->id, 1ul, 0, '\0', return NULL); break;
            case ' ':
                _buf_grow_add_err(kc->id, 1ul, 0, '\0', return NULL);
                if (p == IDTYPE) {
                    char* q = kc->id + h->part[p];
                    if (strncmp(q, "chromosome", strlen(q)) != 0 &&
                            strncmp(q, "supercontig", strlen(q)) != 0 &&
                            strncmp(q, "nonchromosomal", strlen(q)) != 0) {
                        p = UNKNOWN_HDR;
                    }
                } else if (p != ID) {
                    p = UNKNOWN_HDR;
                }
                if (++p < ENS_HDR_PARTCT)
                    h->part[p] = kc->id_l;
                continue;
            case ':':
                _buf_grow_add_err(kc->id, 1ul, 0, '\0', return NULL);
                if (p == ID || p == IDTYPE) {
                    p = UNKNOWN_HDR;
                } else if (p == SEQTYPE) {
                    char* q = kc->id + h->part[p];
                    if (strncmp(q, "dna", strlen(q)) != 0) {
                        if (strncmp(q, "dna_rm", strlen(q)) == 0) {
                            EPR("\nWARNING: reference is repeatmasked");
                        } else if (strncmp(q, "cdna", strlen(q)) == 0) {
                            EPR("\nWARNING: aligning against cDNA");
                        } else if (strncmp(q, "pep", strlen(q)) == 0) {
                            EPR("\nERROR: peptide alignment?");
                            delete h;
                            return NULL;
                        } else if (strncmp(q, "rna", strlen(q)) == 0) {
                            EPR("\nWARNING: aligning against rna");
                        } // else could still be ok.
                    }
                }
                if (++p < ENS_HDR_PARTCT)
                    h->part[p] = kc->id_l;
                continue;
            default:
                if (isspace(c)) {
                    c = '\0';
                    p = UNKNOWN_HDR;
                }
                _buf_grow_add_err(kc->id, 1ul, 0, c, return NULL);
                continue;
        }
        break;
    }
    Mantra contig = {0};
    if (p == NR || p == META) {
        h->p_l = p;
        // ensembl coords are 1-based, we use 0-based.
        contig.corr = atoi(kc->id + h->part[START]) - 1;
        h->end_pos = atoi(kc->id + h->part[END]);
        contig.e = h->end_pos - KEY_WIDTH;
    } else { //TODO: hash lookup from fai
        h->p_l = UNKNOWN_HDR;
        contig.corr = 0;
        contig.e = h->end_pos = ~0u;
        EPQ(dbg > 3, "\nWARNING: non-ensembl reference.");
    }
    h->bnd.push_back(contig);
    return h;
}

static int
end_pos(kct_t*C kc, Hdr* h)
{
    if (h) {
        EPR("processed %lu Nts for %s", kc->s_l - h->s_s + h->bnd.back().corr, kc->id + h->part[0]);

        h->bnd.back().e = kc->s_l - h->s_s + 1;
        h->end_pos = h->bnd.back().e;
        if (dbg > 3)
            show_mantras(kc, h);
    }
    return 0;
}

/*
 * process fasta, store 2bit string and count occurances and locations of keys.
 * store per key the count and the last position.
 */
static int
fa_kc(kct_t* kc, struct gzfh_t* fhin)
{
    void* g;
    int (*gc) (void*);
    uint64_t dna = 0, rc = 0, ndx;
    uint32_t corr = 0;
    int c;
    unsigned i = ~0u; // skip until '>'
    Hdr* h = NULL;
    kc->s_l = 0;
    set_readfunc(fhin, &g, &gc);

    // TODO: realpos based dbsnp or known site boundary insertion.
    while ((c = gc(g)) >= 0) {
        //EPR("%c,%u",c, corr);
        switch(c) {
case 'U': case 'u': c ^= 0x2;
case 'A': case 'a': c ^= 0x3;
case 'T': case 't':
case 'C': case 'c': c ^= 0x2;
case 'G': case 'g': c &= 0x3;
            dna = _seq_next(c, dna, rc);
            _addtoseq(kc->s, c); // kc->s_l grows here.
            if (i == 0) {
                _get_ndx(ndx, ndx, dna, rc);
                if (kc->ndxkct[ndx] != -2u) {
                    ++kc->kct[kc->ndxkct[ndx]];
                } else {
                    _buf_grow(kc->kct, 2, 0);
                    kc->ndxkct[ndx] = kc->kct_l;
                    kc->kct[kc->kct_l++] = 1;
                }
            } else {
                if(--i == KEY_WIDTH - 2) {

                    if (kc->s_l - h->s_s - 1 != 0) { // N-stretch, unless at start, needs insertion
                        end_pos(kc, h);
                        corr += h->bnd.back().corr;
                        h->bnd.push_back({.s = kc->s_l - h->s_s - 1, .e = 0, .corr = 0});
                    }
                    h->bnd.back().corr += corr;
                    corr = 0;
                }
            }
case ' ': case '\t': case '\n': case '\v': case '\f': case '\r':
            break;
case '>':   end_pos(kc, h);
            h = new_header(kc, g, gc);
            ASSERT(h != NULL, return -EFAULT);
            h->s_s = kc->s_l;
            --corr;
default:    ++corr;
            i = KEY_WIDTH - 1;
        }
    }
    ASSERT(h != NULL, return -EFAULT);
    end_pos(kc, h);

    return 0;
}

static void
kct_convert(kct_t* kc)
{
    for (uint64_t ndx = 0ul; ndx != KEYNT_BUFSZ; ++ndx)
        if (kc->ndxkct[ndx] >= kc->kct_l)
            kc->ndxkct[ndx] = kc->kct_l;
}

int
fa_read(struct gzfh_t* fh, kct_t* kc)
{
    int res = -ENOMEM;

    kc->kct = _buf_init_err(kc->kct, 16, goto err);
    kc->id = _buf_init_err(kc->id, 5, goto err);
    kc->s = _buf_init_err(kc->s, 8, goto err);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->ndxkct[i] = -2u;

    /* TODO: load dbSNP and known sites, and mark them. */
    _ACTION(fa_kc(kc, fh + 2), "read and intialized keycounts");
    kct_convert(kc);
    _ACTION(save_seqb2(fh, kc), "writing seqb2: %s", fh[0].name);
    _ACTION(save_kc(fh + 3, kc), "writing keycounts file: %s", fh[3].name);
    _ACTION(reopen(fh + 1, ".nn", ".bd"), "")
    _ACTION(save_boundaries(fh + 1, kc), "writing boundaries: %s", fh[1].name);

    res = 0;
err:
    //fclose occurs in main()
    return res;
}

