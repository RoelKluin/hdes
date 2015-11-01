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
        contig.corr = KEY_WIDTH + atoi(kc->id + h->part[START]) - 1;
        h->end_pos = atoi(kc->id + h->part[END]);
        contig.e = h->end_pos - KEY_WIDTH;
    } else { //TODO: hash lookup from fai
        h->p_l = UNKNOWN_HDR;
        contig.corr = KEY_WIDTH;
        contig.e = h->end_pos = ~0u;
        EPQ(dbg > 3, "\nWARNING: non-ensembl reference.");
    }
    h->bnd.push_back(contig);
    return h;
}

static int
end_pos(kct_t*C kc, Hdr* h, uint32_t corr)
{
    EPR("processed %lu Nts for %s", kc->s_l - h->s_s + corr, kc->id + h->part[0]);
    if (h->end_pos != ~0u) {
        // non-fatal for Y
        EPQ(kc->s_l + corr != h->s_s + h->end_pos, "pos != h->end_pos: %lu == %u",
                kc->s_l - h->s_s + corr, h->end_pos);
        ASSERT(corr != KEY_WIDTH || h->s_s + h->bnd.back().e == kc->s_l, return -EFAULT);
    }
    h->bnd.back().e = kc->s_l - h->s_s + 1;
    h->end_pos = h->bnd.back().e;
    h->end_corr = corr;
    if (dbg > 3)
        show_mantras(kc, h);
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
    int (*ungc) (int, void*);
    uint64_t dna = 0, rc = 0, ndx, b = 0, t = 0ul;
    uint32_t corr = 0u;
    int c, res = -EFAULT;
    Hdr* h = NULL;
    kc->s_l = 0ul;
    set_readfunc(fhin, &g, &gc, &ungc);
    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID
    kc->iter = 0; // here no. boundaries

    // TODO: realpos based dbsnp or known site boundary insertion.
    while (c >= 0) {
        ++kc->iter;
        if (c != '>') { // N-stretch
            do {
                if (c != '\n') {
                    if (c != 'N') {
                        if (B6(b, c) || c == '>') break;
                        EPR("Strange nucleotide:%c, treated as N\n", c);
                    }
                    ++t;
                }
                c = gc(g);
            } while(c != -1);
            // TODO: handle single or a few N's differently.

            if (c == '>' || c == -1) { // skip N's at end.
                corr -= KEY_WIDTH;
                _EVAL(end_pos(kc, h, corr));
                t = res;
		//show_mantras(kc, h);
                continue;
            }
	    if (h->s_s != kc->s_l)
                h->bnd.push_back({0});
            h->bnd.back().s = kc->s_l - h->s_s;
            corr += t;

            h->bnd.back().corr = corr;
        } else {
            // Append the index of the next boundary to the last header.
            // The at_dna is of the last chromo. XXX if at_dna is removed
            // altogether we can stop at h->bnd.end() instead.

            h = new_header(kc, g, gc);
            corr = h->bnd.back().corr;
            h->s_s = kc->s_l;
            ASSERT(h != NULL, return -EFAULT);
            EPQ(dbg > 5, "header %s", kc->id + h->part[0]);
            c = gc(g);
            if (!isb6(b6(c))) {
                t = 0;
                EPR("=>\tN-stretch at Nt %lu", 0);
                continue;
            }
        }
        ungc(c, g);
        for (t = 0; t != KEY_WIDTH; ++t) {
            while (isspace(c = gc(g)));
            b ^= b;
            switch(c) {
case '>':           break;
case 'C': case 'c': b = 2;
case 'G': case 'g': b ^= 1;
case 'U': case 'u':
case 'T': case 't': b ^= 2;
case 'A': case 'a': if (0) {
default:                EPR("FIXME: N's included to make sure the key is completed.\n"
                                "could extend stretch instead?");
                    }
                    dna = _seq_next(b, dna, rc);
                    continue;
            }
            break;
        }
        if_ever (t != KEY_WIDTH) { // c == -1, 'N' or '>'
            // we can get here if an N-stretch follows a new contig
            if (t != 0) {
                // or if two N-s(tretches) lie close to one anorther
                EPQ(dbg > 1, "Key incomplete before next boundary - merging");
            }
            EPQ(dbg > 1, "XXX jump at %lu", kc->s_l - h->s_s + corr);
            t = KEY_WIDTH; // right? [or zero]
            continue; // join regions.
        }
        // ndx (key) is 2nd cNt bit excised, rc or dna
        h->bnd.back().dna = dna;

        while (1) {
            if (isspace(c = gc(g))) continue;
            _get_ndx(ndx, ndx, dna, rc);
            t = kc->ndxkct[ndx];
            if (t != -2u) {
                ++kc->kct[t];
	    } else {
                EPQ(dbg > 5, "new key 0x%lx at %lu", ndx, kc->s_l - h->s_s + corr);
                _buf_grow(kc->kct, 2ul, 0);
                kc->ndxkct[ndx] = kc->kct_l;
		kc->kct[kc->kct_l++] = 1;
            }
            b ^= b;
            switch (c) {
case 'C': case 'c': b = 2;
case 'G': case 'g': b ^= 1;
case 'T': case 't':
case 'U': case 'u': b ^= 2;
case 'A': case 'a':
                // append next 2bit to last ndx. Early Nts are in low bits.
                //EPQ0(dbg > 6, "%s\t%lu\t", kc->id + h->part[0], kc->s_l - h->s_s + corr);
                dna = _seq_next(b, dna, rc);
                _addtoseq(kc->s, b); // kc->s_l grows here.
                continue;
            }
            break;
        }
        if (c == '>' || c == -1) {
            _EVAL(end_pos(kc, h, corr));
            t = res;
        } else {
            // each N-stretch except first adds one extra Nt, for the last key before the Ns.
            EPR("=>\tN-stretch at Nt %lu", kc->s_l - h->s_s + corr);
            t = KEY_WIDTH;
            h->bnd.back().e = kc->s_l - h->s_s;
        }
    }
    _addtoseq(kc->s, 0);
    ASSERT(h != NULL, return -EFAULT);
    res = 0;
err:
    return res;
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

