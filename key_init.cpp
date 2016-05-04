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
#include <unordered_map>
#include <string>
#include <setjmp.h>
#include <cmocka.h> // TODO: unit testing

typedef std::unordered_map<std::string, Hdr*> Hdr_umap;

static inline void
end_pos(kct_t*C kc, Hdr* h)
{
    if (!h) return;
    h->end_pos =  kc->bnd->front().e = kc->s_l - h->s_s;
    kc->totNts += h->end_pos + kc->bnd->front().corr;
    EPR("processed %u(%lu) Nts for %s", h->end_pos + kc->bnd->front().corr, kc->totNts,
            kc->id + h->part[0]);
}

#define ENS_HDR_PARTCT 10
//enum ensembl_parts {ID, SEQTYPE, IDTYPE, IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR};
//ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
// fai does not handle chromosomes with offset.
static Hdr*
new_header(kct_t* kc, Hdr* h, void* g, int (*gc) (void*), Hdr_umap& lookup)
{
    int c;
    uint32_t p = ID, part[ENS_HDR_PARTCT];
    const char* hdr;
    Hdr_umap::const_iterator got;
    part[p] = kc->id_l;

    end_pos(kc, h);

    while ((c = gc(g)) != -1) {
        //EPR("%u\t%c", p, c);
        switch (c) {
            case '\n': _buf_grow_add_err(kc->id, 1ul, 0, '\0', return NULL); break;
            case ' ':
                _buf_grow_add_err(kc->id, 1ul, 0, '\0', return NULL);
                if (p == IDTYPE) {
                    char* q = kc->id + part[p];
                    if (strncmp(q, "chromosome", strlen(q)) != 0 &&
                            strncmp(q, "supercontig", strlen(q)) != 0 &&
                            strncmp(q, "nonchromosomal", strlen(q)) != 0) {
                        p = UNKNOWN_HDR;
                    }
                } else if (p != ID) {
                    p = UNKNOWN_HDR;
                }
                if (p != UNKNOWN_HDR && ++p < ENS_HDR_PARTCT)
                    part[p] = kc->id_l;
                continue;
            case ':':
                _buf_grow_add_err(kc->id, 1ul, 0, '\0', return NULL);
                if (p == ID || p == IDTYPE) {
                    p = UNKNOWN_HDR;
                } else if (p == SEQTYPE) {
                    char* q = kc->id + part[p];
                    if (strncmp(q, "dna", strlen(q)) != 0) {
                        if (strncmp(q, "dna_rm", strlen(q)) == 0) {
                            EPR("\nWARNING: reference is repeatmasked");
                        } else if (strncmp(q, "cdna", strlen(q)) == 0) {
                            EPR("\nWARNING: aligning against cDNA");
                        } else if (strncmp(q, "pep", strlen(q)) == 0) {
                            EPR("\nERROR: peptide alignment?");
                            return NULL;
                        } else if (strncmp(q, "rna", strlen(q)) == 0) {
                            EPR("\nWARNING: aligning against rna");
                        } // else could still be ok.
                    }
                }
                if (p != UNKNOWN_HDR && ++p < ENS_HDR_PARTCT)
                    part[p] = kc->id_l;
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
    if (c == -1)
        return NULL;
    hdr = kc->id + part[ID];
    got = lookup.find(hdr);

    if (got == lookup.end()) {

        //    new Hdr;
        _buf_grow_err(kc->h, 1ul, 0, return NULL);
        h = kc->h + kc->h_l++;

        h->part = (uint32_t*)malloc(++p * sizeof(uint32_t));
        NB(h->part != NULL);
        memcpy(h->part, part, p * sizeof(*part));

        std::pair<std::string,Hdr*> hdr_entry(hdr, h);
        lookup.insert(hdr_entry);
        h->s_s = kc->s_l;
        kc->bnd->push_front({0});
    } else {

        EPR("contig occurred twice: %s", hdr);
        // To fix order reversal would need kc->s movement and adaptations including h->s_s.
        NB(h->part[START] > got->second->part[END], "Duplicate entries in reversed order");
        h = got->second;
        NB(h == kc->h + kc->h_l - 1, "Duplicate entries, but not in series");

        // only insert when last contig had any sequence
        if (kc->bnd->front().e != kc->bnd->front().s) {
            if (p != NR && p != META) {
                h->p_l = UNKNOWN_HDR;
                h->end_pos = ~0u;
                WARN("No offsets recognized in 2nd header, sequence will be concatenated.");
                return h;
            }

            // correction propagation Not thoroughly checked yet...
            pos_t corr = kc->bnd->front().corr - kc->bnd->front().e; // start of current is added later.
            kc->bnd->push_front({.s=0, .e=0, .corr=corr });
        }
    }

    if (p == NR || p == META) {
        h->p_l = p;
        // ensembl coords are 1-based, we use 0-based.
        kc->bnd->front().corr += atoi(kc->id + h->part[START]) - 1;
        h->end_pos = atoi(kc->id + h->part[END]);
        kc->bnd->front().e = h->end_pos;
    } else { //TODO: hash lookup from fai
        h->p_l = UNKNOWN_HDR;
        kc->bnd->front().corr = 0;
        kc->bnd->front().e = h->end_pos = ~0u;
        EPR("\nWARNING: non-ensembl reference.");
    }
    return h;
}

/*
 * process fasta, store 2bit string and count occurances and locations of keys.
 * store per key whether it occurred multiple times the last position.
 */
static int
fa_kc(kct_t* kc, struct gzfh_t* fhin)
{
    void* g;
    int (*gc) (void*);
    pos_t corr = 0;
    unsigned i = ~0u; // skip until '>'
    Hdr* h = NULL;
    keyseq_t seq = {.p = 0, .dna = 0x000000ff, .rc = 0x000000ff, .t = 0x0000ffff};
    kc->s_l = 0;
    kc->uqct = 0;
    kc->totNts = 0;
    kc->bnd = new std::forward_list<Mantra>();
    Hdr_umap lookup;
    HK hk = {.hoffs = kc->h_l};
    set_readfunc(fhin, &g, &gc);

    // TODO: realpos based dbsnp or known site boundary insertion.
    while ((seq.t = gc(g)) <= 0xff) {
        switch((seq.t | i) & ~0x20) {
case 'U':   seq.t ^= 0x2;
case 'A':   seq.t ^= 0x3;
case 'T':
case 'C':   seq.t ^= 0x2;
case 'G':   seq.t &= 0x3;
        {
            seq_next(seq);
            seq_t ndx;
            _addtoseq(kc->s, seq.t); // kc->s_l grows here.
            seq_t* n = kc->contxt_idx + get_kct0(kc, seq, ndx);
            if (*n == NO_KCT) {
                _buf_grow(kc->kct, 2, 0);
                *n = kc->kct_l++;
                // set first pos + orient
                kc->kct[*n] = (kc->s_l - h->s_s) << 1 | (seq.t != 0);
                NB(b2pos_of(kc->kct[*n]) > KEY_WIDTH - 1);
            } else {
                if (!(kc->kct[*n] & DUP_BIT)) {
                    kc->kct[*n] |= DUP_BIT;   // mark it as dup
                    --kc->uqct;
                }
            }
            break;
        }
case 'N':   i = (KEY_WIDTH - 1) << 8;
            //EPR("%seq.t",seq.t);
case 'N' | ((KEY_WIDTH - 1) << 8):
            ++corr;
            break;
default:    if (isspace(seq.t))
                continue;
            switch (seq.t & ~0x20) {
    case 'U':   seq.t ^= 0x2;
    case 'A':   seq.t ^= 0x3;
    case 'T':
    case 'C':   seq.t ^= 0x2;
    case 'G':   seq.t &= 0x3;
                EPR0(".");
                if (i == ((KEY_WIDTH - 1) << 8)) { // key after header/stretch to be rebuilt
                    NB(h != NULL);
                    if (kc->s_l != h->s_s) { // N-stretch, unless at start, needs insertion
                        end_pos(kc, h);
                        corr += kc->bnd->front().corr;
                        NB(h->s_s < kc->s_l);
                        NB(kc->s_l - h->s_s <= 0x3fffffff,"TODO: split seq for huge contigs");
                        kc->bnd->push_front({.s = (pos_t)(kc->s_l - h->s_s)});
                    }
                    kc->bnd->front().corr += corr;
                    corr = 0;
                }
                i -= 0x100;
                seq_next(seq);
                _addtoseq(kc->s, seq.t);
                //print_dna(seq.dna);
                break;
    case 0x1e:{ if (h) {
                    hk.koffs = kc->kct_l;
                    _buf_grow_add_err(kc->hk, 1ul, 0, hk, return -ENOMEM);
                }
                h = new_header(kc, h, g, gc, lookup);
                NB(h != NULL);

                i = (KEY_WIDTH - 1) << 8;
                corr = -1u;
            }
    default:    ++corr;
            }
        }
    }
    kc->uqct += kc->kct_l;
    NB(h != NULL);
    hk.koffs = kc->kct_l;
    _buf_grow_add_err(kc->hk, 1ul, 0, hk, return -ENOMEM);
    end_pos(kc, h);
    fprintf(stderr, "Initial unique keys: %u / %u\n", kc->uqct, kc->kct_l);
    return 0;
}

int
fa_read(struct gzfh_t* fh, kct_t* kc)
{
    int res = -ENOMEM;

    kc->kct = _buf_init_err(kc->kct, 16, goto err);
    kc->id = _buf_init_err(kc->id, 5, goto err);
    kc->s = _buf_init_err(kc->s, 8, goto err);
    kc->h = _buf_init_err(kc->h, 1, goto err);
    kc->hk = _buf_init_err(kc->hk, 1, goto err);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->contxt_idx[i] = NO_KCT;

    /* TODO: load dbSNP and known sites, and mark them. */
    _ACTION(fa_kc(kc, fh + 2), "read and intialized keycounts");
    _ACTION(save_seqb2(fh, kc), "writing seqb2: %s", fh[0].name);
    _ACTION(save_kc(fh + 3, kc), "writing keycounts file: %s", fh[3].name);
    _ACTION(reopen(fh + 1, ".nn", ".bd"), "");
    _ACTION(save_boundaries(fh + 1, kc), "writing boundaries: %s", fh[1].name);

    res = 0;
err:
    //fclose occurs in main()
    return res;
}




