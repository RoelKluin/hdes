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
parse_header_parts(Key_t* kc, gzin_t* gz, uint32_t* part)
{
    int p = ID;
    part[p] = kc->id_l;
    int c;
    while ((c = GZTC(gz)) != -1) {
        //EPR("%u\t%c", p, c);
        switch (c) {
            case '\n': buf_grow_add(kc->id, 1, '\0'); break;
            case ' ':
                buf_grow_add(kc->id, 1, '\0');
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
                buf_grow_add(kc->id, 1, '\0');
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
                            return -EFAULT;
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
                buf_grow_add(kc->id, 1, c);
                continue;
        }
        break;
    }
    return p == NR || p == META ? p : UNKNOWN_HDR;
}

typedef std::unordered_map<std::string, Hdr*> Hdr_umap;

//enum ensembl_part {ID, SEQTYPE, IDTYPE, IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR};
//ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
// fai does not handle chromosomes with offset.
static Hdr*
new_header(Key_t* kc, Hdr* h, gzin_t* gz, Hdr_umap& lookup)
{
    uint32_t part[ENS_HDR_PARTCT] = {0};
    int res = parse_header_parts(kc, gz, part);
    if (res < 0)
        return NULL;

    const char* hdr = kc->id + part[ID];
    Hdr_umap::const_iterator got = lookup.find(hdr);
    unsigned offs = 0;

    if (got == lookup.end()) {

        //    new Hdr;
        buf_grow_struct_add(kc->h, 0); // Note: initializes all members to zero or NULL
	h = kc->h + (kc->h_l - 1);
	h->ido = part[ID];

        std::pair<std::string,Hdr*> hdr_entry(hdr, h);
        lookup.insert(hdr_entry);

        buf_grow_struct_add(kc->bnd, .ho = (uint32_t)(h - kc->h), .s = NT_WIDTH);
	// ensembl coords are 1-based, we use 0-based.
	if (res == NR || res == META)
	    offs =  atoi(kc->id + part[START]) - 1;
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
        offs = kc->bnd[kc->bnd_l-1].s;
        if (offs != 0) {
            if (res == NR || res == META) {

                buf_grow_struct_add(kc->bnd, .ho = (uint32_t)(h - kc->h), .s= NT_WIDTH);
		offs = h->noffs[h->noffs_l - 1].corr;
		res = h->p_l; // keep status observed in 1st header
	    } else {
	        WARN("No offsets recognized in 2nd header, sequence will be concatenated.");
	    }
        }
    }
    h->p_l = res;
    buf_grow_struct_add(h->noffs, .pos = 0, .corr = offs);
    return h;
}

static inline int
finish_contig(Key_t*C kc, Hdr* h, uint32_t p)
{
    // the 2bit buffer per contig starts at the first nt 0 of 4.
    h->len = (p >> 3) + !!(p & 6);
    h->end = (p >> 1) + h->noffs[h->noffs_l-1].corr;

    buf_grow_add(kc->hkoffs, 1, kc->kct_l);
    kc->bnd[kc->bnd_l-1].e = p + 2;
    kc->totNts += h->end;
    EPR("Processed %u Nts for %s => %lu", h->end, kc->id + h->ido, kc->totNts);
    if (kc->s_l & 3) {
        // the 2bit buffer per contig starts at the first nt 0 of 4.
        kc->s_l += -kc->s_l & 3;
        buf_grow_shift(kc->s, 1ul, 2);
        kc->s[kc->s_l>>2] = '\0';
    }
    return 0;
}

static inline void
addtoseq(Key_t* kc, keyseq_t& seq)
{
    seq.rc = ((seq.rc << 2) & KEYNT_MASK) | (seq.t ^ 2);
    seq.dna = seq.t << KEYNT_TOP | seq.dna >> 2;
    seq.p += 2; // first bit is for orientation
    if (kc->s_l & 3) {
        seq.t = kc->s[kc->s_l>>2] | (seq.t << ((kc->s_l & 3) << 1));
    } else {
        buf_grow_shift(kc->s, 1ul, 2);
    }
    kc->s[kc->s_l++>>2] = seq.t;
}

/*
 * process fasta, store 2bit string and count occurances and locations of keys.
 * store per key whether it occurred multiple times the last position.
 */
static int
fa_kc(Key_t* kc, gzin_t* gz)
{
    uint32_t N_count = 0;
    unsigned i = ~0u; // skip until '>'
    int res;
    Hdr* h = NULL;
    keyseq_t seq = {0};
    Hdr_umap lookup;

    // TODO: realpos based dbsnp or known site boundary insertion.
    while ((seq.t = GZTC(gz)) <= 0xff) {
        switch(seq.t & ~0x20) {
case 'U':   seq.t ^= 0x2;
case 'A':   seq.t ^= 0x3;
case 'T':
case 'C':   seq.t ^= 0x2;
case 'G':   seq.t &= 0x3;
            addtoseq(kc, seq);
            if (i == 0) {
                uint32_t* n = get_kct(kc, seq, 1);
                if (*n == NO_K) {

                    buf_grow(kc->kct, 1);
                    *n = kc->kct_l++;
                    // set first pos + orient, one based & one left shifted pos
                    kc->kct[*n] = seq.p;

                } else if (!(kc->kct[*n] & DUP_BIT)) {

                    --kc->ct;
                    kc->kct[*n] |= DUP_BIT;   // mark it as dup
                }
                seq.p = b2pos_of(seq.p);
            } else {
                if (i-- == KEY_WIDTH - 1) { // key after header/stretch to be rebuilt
                    NB(h != NULL);
                    if (seq.p > 2u) { // N-stretch, unless at start, needs insertion
                        NB((seq.p & 1) == 0);
                        kc->bnd[kc->bnd_l-1].e = seq.p;
			uint32_t prev_corr = h->noffs[h->noffs_l-1].corr;
                        kc->totNts += (seq.p >> 1) + prev_corr;

			uint32_t pos = seq.p + NT_WIDTH - 2;
                        buf_grow_struct_add(h->noffs, .pos = pos, .corr = N_count + prev_corr);
			buf_grow_struct_add(kc->bnd, .ho = (uint32_t)(h - kc->h), .s = pos);
                    }
                    N_count = 0;
                }
            }
            break;
case 0x1e: // new contig
            NB(seq.t == '>');
            if (h) {
                // FIXME: Y-contig occurs twice. Need to lookup here and skip if already present.
                _EVAL(finish_contig(kc, h, seq.p));
                seq.p = 0;
            }
            h = new_header(kc, h, gz, lookup);
            NB(h != NULL);
            i = KEY_WIDTH - 1;
            break;
default:    if (isspace(seq.t))
                break;
            i = KEY_WIDTH - 1;
            ++N_count;
        }
    }
    res = finish_contig(kc, h, seq.p);
    kc->ct += kc->kct_l;
    fprintf(stderr, "Initial unique keys: %u / %u\n", kc->ct, kc->kct_l);
err:
    return res;
}

int
fa_read(struct gzfh_t* fh, Key_t* kc)
{
    int res = -ENOMEM;
    gzin_t gz = {0};

    kc->kct = buf_init(kc->kct, 8);
    kc->id = buf_init(kc->id, 8);
    kc->s = buf_init(kc->s, 8);
    kc->h = buf_init(kc->h, 1);
    kc->hkoffs = buf_init(kc->hkoffs, 1);
    kc->bnd = buf_init(kc->bnd, 8);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->contxt_idx[i] = NO_K;

    /* TODO: load dbSNP and known sites, and mark them. */
    set_readfunc(fh + 2, &gz);
    _ACTION(fa_kc(kc, &gz), "read and intialized keycounts");
    _ACTION(save_seqb2(fh + 1, kc), "writing seqb2: %s", fh[1].name);
    _ACTION(save_kc(fh, kc), "writing keycounts file: %s", fh[0].name);

    res = 0;
err:
    //fclose occurs in main()
    return res;
}




