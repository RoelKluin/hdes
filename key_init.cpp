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

char*
get_header_part(char *s, ensembl_part tgt)
{
    NB(tgt != ID);
    for (unsigned i = ID; *s != '\0' || ++i != tgt; ++s)
        {}
    return s;
}

#define ENS_HDR_PARTCT 10
static int
parse_header_parts(kct_t* kc, void* g, int (*gc) (void*), uint32_t* part)
{
    int p = ID;
    part[p] = kc->id_l;
    int c;
    while ((c = gc(g)) != -1) {
        //EPR("%u\t%c", p, c);
        switch (c) {
            case '\n': buf_grow_add(kc->id, 1ul, 0, '\0'); break;
            case ' ':
                buf_grow_add(kc->id, 1ul, 0, '\0');
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
                buf_grow_add(kc->id, 1ul, 0, '\0');
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
                buf_grow_add(kc->id, 1ul, 0, c);
                continue;
        }
        break;
    }
    return p;
}

static inline void
set_header_type(kct_t*C kc, Hdr* h, int type, uint32_t corr)
{
    h->p_l = type;
    kc->bnd[kc->bnd_l-1].corr = corr;
}

typedef std::unordered_map<std::string, Hdr*> Hdr_umap;

//enum ensembl_part {ID, SEQTYPE, IDTYPE, IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR};
//ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
// fai does not handle chromosomes with offset.
static Hdr*
new_header(kct_t* kc, Hdr* h, void* g, int (*gc) (void*), Hdr_umap& lookup, uint32_t endpos)
{
    uint32_t part[ENS_HDR_PARTCT] = {0};
    int res = parse_header_parts(kc, g, gc, part);
    if (res < 0)
        return NULL;

    const char* hdr = kc->id + part[ID];
    Hdr_umap::const_iterator got = lookup.find(hdr);

    if (got == lookup.end()) {

        //    new Hdr;
        buf_grow(kc->h, 1ul, 0);
        h = kc->h + kc->h_l++;
        h->ido = part[ID];

        std::pair<std::string,Hdr*> hdr_entry(hdr, h);
        lookup.insert(hdr_entry);

        buf_grow(kc->bnd, 1ul, 0);
        kc->bnd[kc->bnd_l++] = {.ho = h - kc->h, .s= NT_WIDTH};
    } else {

        EPR("contig occurred twice: %s", hdr);
        // To fix order reversal would need kc->s movement and adaptations in header start.
        if (h->p_l != UNKNOWN_HDR && got->second->p_l != UNKNOWN_HDR) {

            NB(atoi(kc->id + part[START]) > atoi(get_header_part(kc->id + got->second->ido, END)),
                    "Duplicate entries in reversed order");
        }
        h = got->second;
        NB(h == kc->h + kc->h_l - 1, "Duplicate entries, but not in series");

        // only insert when last contig had at least one KEY_WIDTH of sequence
        if (endpos != kc->bnd[kc->bnd_l-1].s) {
            if (res != NR && res != META) {
                set_header_type(kc, h, UNKNOWN_HDR, kc->bnd[kc->bnd_l-1].corr);
                WARN("No offsets recognized in 2nd header, sequence will be concatenated.");
                return h;
            }

            // correction propagation Not thoroughly checked yet...
            uint32_t corr = kc->bnd[kc->bnd_l-1].corr - endpos; // start of current is added later.
            buf_grow(kc->bnd, 1ul, 0);
            kc->bnd[kc->bnd_l++] = {.ho = h - kc->h, .s= NT_WIDTH, .corr=corr};
        }
    }

    if (res == NR || res == META) {
        // ensembl coords are 1-based, we use 0-based.
        set_header_type(kc, h, res, atoi(kc->id + part[START]) - 1);
    } else { //TODO: hash lookup from fai
        set_header_type(kc, h, UNKNOWN_HDR, 0);
    }
    return h;
}

static inline void
end_pos(kct_t*C kc, uint32_t len)
{
    kc->bnd[kc->bnd_l-1].e = len + 2;
    kc->totNts += len + kc->bnd[kc->bnd_l-1].corr;
}

static inline int
finish_contig(kct_t*C kc, Hdr* h, keyseq_t &seq)
{
    // the 2bit buffer per contig starts at the first nt 0 of 4.
    h->len = (seq.p >> 3) + !!(seq.p & 6);
    h->end = seq.p;
    buf_grow_add(kc->hkoffs, 1ul, 0, kc->kct_l);
    end_pos(kc, seq.p);
    EPR("processed %u(%lu) Nts for %s", seq.p >> 1, kc->totNts, kc->id + h->ido);
    return 0;
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
    uint32_t corr = 0;
    unsigned i = ~0u; // skip until '>'
    int res;
    Hdr* h = NULL;
    keyseq_t seq = {0};
    Hdr_umap lookup;
    set_readfunc(fhin, &g, &gc);

    // TODO: realpos based dbsnp or known site boundary insertion.
    while ((seq.t = gc(g)) <= 0xff) {
        switch(seq.t & ~0x20) {
case 'U':   seq.t ^= 0x2;
case 'A':   seq.t ^= 0x3;
case 'T':
case 'C':   seq.t ^= 0x2;
case 'G':   seq.t &= 0x3;
            _addtoseq(kc->s, seq);
            next_seqpos(kc->s, seq);
            if (i == 0) {
                uint32_t* n = get_kct(kc, seq, 1);
                if (*n == NO_K) {

                    buf_grow(kc->kct, 1, 0);
                    *n = kc->kct_l++;
                    // set first pos + orient, one based & one left shifted pos
                    kc->kct[*n] = seq.p;

                } else if (!(kc->kct[*n] & DUP_BIT)) {

                    kc->kct[*n] |= DUP_BIT;   // mark it as dup
                    --kc->ct;
                }
                seq.p = b2pos_of(seq.p);
            } else {
                if (i-- == KEY_WIDTH - 1) { // key after header/stretch to be rebuilt
                    NB(h != NULL);
                    if (seq.p > 2u) { // N-stretch, unless at start, needs insertion
                        NB((seq.p & 1) == 0);
                        end_pos(kc, seq.p - 2);
                        corr += kc->bnd[kc->bnd_l-1].corr;
                        buf_grow(kc->bnd, 1ul, 0);
                        kc->bnd[kc->bnd_l++] = {.ho = h - kc->h, .s = seq.p + NT_WIDTH - 2};
                    }
                    kc->bnd[kc->bnd_l-1].corr += corr;
                    corr = 0;
                }
            }
            break;
case 0x1e: // new contig
            NB(seq.t == '>');
            if (h) {
                // FIXME: Y-contig occurs twice. Need to lookup here and skip if already present.
                _EVAL(finish_contig(kc, h, seq));
                if (kc->s_l & 3) {
                    // the 2bit buffer per contig starts at the first nt 0 of 4.
                    kc->s_l += -kc->s_l & 3;
                    buf_grow(kc->s, 1ul, 2);
                    kc->s[kc->s_l>>2] = '\0';
                }
                seq.p = 0;
            }
            h = new_header(kc, h, g, gc, lookup, seq.p);
            NB(h != NULL);
            i = KEY_WIDTH - 1;
            break;
default:    if (isspace(seq.t))
                break;
            i = KEY_WIDTH - 1;
            ++corr;
        }
    }
    res = finish_contig(kc, h, seq);
    kc->ct += kc->kct_l;
    fprintf(stderr, "Initial unique keys: %u / %u\n", kc->ct, kc->kct_l);
err:
    return res;
}

int
fa_read(struct gzfh_t* fh, kct_t* kc)
{
    int res = -ENOMEM;

    kc->kct = buf_init(kc->kct, 8);
    kc->id = buf_init(kc->id, 8);
    kc->s = buf_init(kc->s, 8);
    kc->h = buf_init(kc->h, 1);
    kc->hkoffs = buf_init(kc->hkoffs, 1);
    kc->bnd = buf_init(kc->bnd, 8);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->contxt_idx[i] = NO_K;

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




