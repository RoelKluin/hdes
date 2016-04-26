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

typedef std::unordered_map<std::string, Hdr*> Hdr_umap;

static inline void
end_pos(kct_t*C kc, Hdr* h)
{
    if (h) {
        h->end_pos =  h->bnd->back().e = kc->s_l - h->s_s;
        kc->totNts += h->end_pos + h->bnd->back().corr;
        EPR("processed %u(%lu) Nts for %s", h->end_pos + h->bnd->back().corr, kc->totNts,
                kc->id + h->part[0]);

        if (dbg > 3)
            show_mantras(kc, h);
    }
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
            case '\n': _buf_grow_add_err(kc->id, 1ul, 0, '\0', h = NULL; goto err); break;
            case ' ':
                _buf_grow_add_err(kc->id, 1ul, 0, '\0', h = NULL; goto err);
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
                _buf_grow_add_err(kc->id, 1ul, 0, '\0', h = NULL; goto err);
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
                _buf_grow_add_err(kc->id, 1ul, 0, c, h = NULL; goto err);
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
        ASSERT(h->part != NULL, return NULL);
        memcpy(h->part, part, p * sizeof(*part));

        h->bnd = new std::list<Mantra>();
        std::pair<std::string,Hdr*> hdr_entry(hdr, h);
        lookup.insert(hdr_entry);
        h->s_s = kc->s_l;
        h->bnd->push_back({0});
    } else {

        EPR("contig occurred twice: %s", hdr);
        // To fix order reversal would need kc->s movement and adaptations including h->s_s.
        ASSERT(h->part[START] > got->second->part[END], h = NULL; goto err, "Duplicate entries in reversed order");
        h = got->second;
        ASSERT(h == kc->h + kc->h_l - 1, h = NULL; goto err, "Duplicate entries, but not in series");

        // only insert when last contig had any sequence
        if (h->bnd->back().e != h->bnd->back().s) {
            if (p != NR && p != META) {
                h->p_l = UNKNOWN_HDR;
                h->end_pos = ~0u;
                WARN("No offsets recognized in 2nd header, sequence will be concatenated.");
                return h;
            }

            // correction propagation Not thoroughly checked yet...
            pos_t corr = h->bnd->back().corr - h->bnd->back().e; // start of current is added later.
            h->bnd->push_back({.s=0, .e=0, .corr=corr });
        }
    }

    if (p == NR || p == META) {
        h->p_l = p;
        // ensembl coords are 1-based, we use 0-based.
        h->bnd->back().corr += atoi(kc->id + h->part[START]) - 1;
        h->end_pos = atoi(kc->id + h->part[END]);
        h->bnd->back().e = h->end_pos - KEY_WIDTH;
    } else { //TODO: hash lookup from fai
        h->p_l = UNKNOWN_HDR;
        h->bnd->back().corr = 0;
        h->bnd->back().e = h->end_pos = ~0u;
        EPQ(dbg > 3, "\nWARNING: non-ensembl reference.");
    }
err:
    if (h == NULL)
        delete h->bnd;
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
    int c;
    unsigned i = ~0u; // skip until '>'
    seq_t ndx;
    Hdr* h = NULL;
    keyseq_t seq = {0};
    HK hk = {0};
    hk.koffs;
    kc->s_l = 0;
    kc->uqct = 0;
    kc->totNts = 0;
    Hdr_umap lookup;
    set_readfunc(fhin, &g, &gc);

    // TODO: realpos based dbsnp or known site boundary insertion.
    while ((c = gc(g)) >= 0) {
        //EPR("%c," Pfmt,c, corr);
        switch((c | i) & ~0x20) {
case 'U':   c ^= 0x2;
case 'A':   c ^= 0x3;
case 'T': 
case 'C':   c ^= 0x2;
case 'G':   c &= 0x3;
        {
            seq.dna = _seq_next(c, seq);
            //print_dna(seq.dna);
            _addtoseq(kc->s, c); // kc->s_l grows here.
            seq_t* n = kc->ndxkct + _get_kct0(kc, seq, seq.t, ndx, return -EFAULT);
            if (*n == NO_KCT) {
                ++hk.koffs;
                _buf_grow(kc->kct, 2, 0);
                *n = kc->kct_l++;
                kc->kct[*n] = 0ul;
                ++kc->uqct;
            } else {
                C uint64_t m = 1ul << ORIENT_SHFT;
                kc->kct[*n] &= m ^ -m; // clear position and orientation
                if (!(kc->kct[*n] & DUP_BIT)) {
                    kc->kct[*n] |= DUP_BIT;                    // mark it's a dup
                    --kc->uqct;
                }
            }
            EPQ(dbg >4, "[%lu, " Sfmt ", 0x%lx, " Sfmt "]:debug %u", kc->s_l, seq.dna,
                    n - kc->ndxkct, *n, print_dna(seq.dna));
            /*ASSERT((kc->kct[*n] & ((seq.t != 0) << ORIENT_SHFT)) == 0 &&
                    (kc->kct[*n] & (kc->s_l + 1)) == 0 &&
                    (((seq.t != 0) << ORIENT_SHFT) & (kc->s_l + 1)) == 0, return -EFAULT);*/
            kc->kct[*n] ^= ((uint64_t)(seq.t != 0) << ORIENT_SHFT) | (kc->s_l - h->s_s); // set latest pos + orient
            break;
        }
case 'N':   i = (KEY_WIDTH - 1) << 8;
            //EPR("%c",c);
case 'N' | ((KEY_WIDTH - 1) << 8):
            ++corr;
            break;
default:    if (isspace(c))
                continue;
            switch (c & ~0x20) {
    case 'U':   c ^= 0x2;
    case 'A':   c ^= 0x3;
    case 'T': 
    case 'C':   c ^= 0x2;
    case 'G':   c &= 0x3;
                if (i == ((KEY_WIDTH - 1) << 8)) {
                    ASSERT(h, return -EFAULT);
                    if (kc->s_l != h->s_s) { // N-stretch, unless at start, needs insertion
                        end_pos(kc, h);
                        corr += h->bnd->back().corr;
                        ASSERT(h->s_s > kc->s_l, return -EFAULT);
                        ASSERT(kc->s_l - h->s_s <= 0xffffffff, return -EFAULT);
                        h->bnd->push_back({.s = (pos_t)(kc->s_l - h->s_s), .e = 0, .corr = 0});
                    }
                    h->bnd->back().corr += corr;
                    corr = 0;
                }
                i -= 0x100;
                seq.dna = _seq_next(c, seq);
                _addtoseq(kc->s, c);
                //print_dna(seq.dna);
                break;
    case 0x1e:  if (h) {
                    // FIXME: convert header to reallocable for this purpose
                    hk.hoffs = kc->h_l - 1;
                    _buf_grow_add_err(kc->hk, 1ul, 0, hk, return -ENOMEM);
                }
                h = new_header(kc, h, g, gc, lookup);
                // hk.koffs is cumulutive;
                ASSERT(h != NULL, return -EFAULT);

                i = (KEY_WIDTH - 1) << 8;
                corr = -1u;
    default:    ++corr;
            }
        }
    }
    ASSERT(h != NULL, return -EFAULT);
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
        kc->ndxkct[i] = NO_KCT;

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

