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
    h->bnd.back().e = kc->s_l - h->s_s;
    h->end_pos = h->bnd.back().e;
    h->end_corr = corr;
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
            h->bnd.back().s = kc->s_l - h->s_s;
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

            corr += t;
            if (c == '>' || c == -1) { // skip N's at end.
                corr -= KEY_WIDTH;
                h->bnd.back().corr = corr;
                _EVAL(end_pos(kc, h, corr));
                t = res;
                continue;
            }

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
default: // include 'N's and such to make sure the key is completed. FIXME: extend stretch instead?
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
            uint64_t *ct;
            _get_ndx(ndx, ndx, dna, rc);
            if (kc->ndxkct[ndx] != -2u) {
                t = kc->ndxkct[ndx];
                if (ndx == dbgndx) dbgkct = t;
                //dbg = t == dbgkct ? dbg | 8 : dbg & ~8;
                ASSERT(t < kc->kct_l, return -EFAULT, "0x%lx, 0x%lx", t, ndx);
                ct = kc->kct + t + 1;
                if ((*ct & B2LEN_MASK) != B2LEN_MASK) {
                    EPQ(dbg > 7, "sequence & length");
                    t = *ct >> B2LEN_SHFT;

                    *ct += ONE_B2SEQ;
                    if (t < 32ul) {
                        --ct;
                    } else if (t == 61ul) {
                        EPQ(dbg > 7, "conversion to index");
                        t = *ct & B2SEQ_MASK;  // temp store 2bit dna
                        *ct = B2LEN_MASK | 62; // conversion: set to index, make sure B2LEN_MASK bits are set
                        b = *--ct;
                        *ct = (uint64_t)malloc(2 * sizeof(uint64_t)); // for 2bits
                        ct = (uint64_t*)*ct;
                        *ct = b;                // added after current 2bit dna, in ndx
                        *++ct = t;
                        t = 61;
                    }
                } else {
                    EPQ(dbg > 7, "in index format:");
                    t = (*ct--)++; // incr length one, move to address.
                    ct = (uint64_t *)*ct;

                    if ((t & 0x3f) == 0ul) {
                        if ((t & 0xffffffff) == (0x40u << ((t >> 32) & 0xff))) {
                            t = kc->kct[kc->ndxkct[ndx] + 1] += 1ul << 32; // increment m
                            ct = (uint64_t *)realloc(ct,
                                    (2 * sizeof(uint64_t)) << ((t >> 32) & 0xff));
                            ASSERT(ct != NULL, return -ENOMEM);
                            kc->kct[kc->ndxkct[ndx]] = (uint64_t)ct;
                            --t;
                        }
                        ct += (t & 0xffffffff) >> 5; // point to the u64 to be updated
                        *ct = ct[1] = 0ul;
                    } else {
                        ct += (t & 0xffffffff) >> 5; // point to the u64 to be updated
                    }
                }
            } else {
                EPQ(dbg > 5, "new key 0x%lx at %lu", ndx, kc->s_l - h->s_s + corr);
                _buf_grow(kc->kct, 2ul, 0);
                kc->ndxkct[ndx] = kc->kct_l;
                ct = kc->kct + kc->kct_l;
                kc->kct_l += 2;
                *ct = t = 0ul;
                ct[1] = ONE_B2SEQ;
            }
            b ^= b;
            switch (c) {
case 'C': case 'c': b = 2;
case 'G': case 'g': b ^= 1;
case 'T': case 't':
case 'U': case 'u': b ^= 2;
case 'A': case 'a':
                // append next 2bit to last ndx. Early Nts are in low bits.
                dna = _seq_next(b, dna, rc);
                *ct |= b << ((t & 0x1f) << 1);
                EPQ0(dbg > 6, "%s\t%lu\t", kc->id + h->part[0], kc->s_l - h->s_s + corr);
                print_dna(*ct, dbg > 6, '\n', (ct >= kc->kct) && (ct < (kc->kct + kc->kct_l)) ? 29 : 32);
                _addtoseq(kc->s, b); // kc->s_l grows here.
                continue;
            }
            /*if ((ct >= kc->kct) && (ct < (kc->kct + kc->kct_l))) { // pointer in range - valid for non-void pointers
                if (t < 32) ++ct;
                *ct -= ONE_B2SEQ;
            } else {// If not within kct range then this is an extended keycount:
                --kc->kct[kc->ndxkct[ndx] + 1];
            }*/
            // at a boundary no next Nt is inserted, but its position is updated. it
            // won't be moved.
            _addtoseq(kc->s, 0);
            break;
        }

        h->bnd.back().end_dna = dna;
        if (c == '>' || c == -1) {
            _EVAL(end_pos(kc, h, corr));
            t = res;
        } else {
            // each N-stretch except first adds one extra Nt, for the last key before the Ns.
            EPR("=>\tN-stretch at Nt %lu", kc->s_l - h->s_s + corr);
            t = KEY_WIDTH;
            h->bnd.back().e = kc->s_l - h->s_s;
            h->bnd.push_back({0});
        }
    }
    //_addtoseq(kc->s, 0);
    ASSERT(h != NULL, return -EFAULT);
    res = 0;
err:
    return res;
}

// concatenate 2bits, byte aligned, for each key in buffer kc->ts,
// store lengths in kc->kct;
static int
kct_convert(kct_t* kc)
{
    kc->ts_l = kc->s_l; // about as many next NTs as Nts.
    // XXX: requirement of addition below is strangely needed +1 seems to work in most cases. +4 here just to be sure..
    // Maybe because end of contig key also needs an nextnt to keep its counts? in that case the
    // addition should be the no. contigs + N-stretches ..?

    const uint64_t ts_end = (kc->ts_l >> 2) + !!(kc->ts_l & 3) + kc->iter;
    uint8_t *dest = (uint8_t*)malloc(ts_end);
    ASSERT(dest != NULL, return -ENOMEM);
    kc->ts = dest;
    *dest = '\0';
    uint64_t offs = 0ul;

    for (uint64_t ndx = 0ul; ndx != KEYNT_BUFSZ; ++ndx) {
        if (kc->ndxkct[ndx] >= kc->kct_l) {
            kc->ndxkct[ndx] = kc->kct_l;
            continue;
        }
        uint64_t *src = kc->kct + kc->ndxkct[ndx];
        dbg = kc->ndxkct[ndx] == dbgkct ? dbg | 8 : dbg & ~8;
        uint64_t l = src[1];
        unsigned __int128 x;
        unsigned t = offs & 3;
        ASSERT(dest < kc->ts + ts_end, return -EFAULT);
        if ((l & B2LEN_MASK) != B2LEN_MASK) { // then x contains count and next Nts below this bit
            x = ((unsigned __int128)(l & B2SEQ_MASK) << 64) | src[0];

            l >>= B2LEN_SHFT;            // isolate twobit count
            EPQ(dbg > 3 && l == 0, "No nextNts for 0x%x", kc->ndxkct[ndx]); // can happen at boundaries
            src[1] = (l << BIG_SHFT) | offs;
            dest = kc->ts + (offs >> 2); // destination for copy.
            offs += l;                   // update offs for next round.
            x <<= t << 1;                // prealign 2bits for appending to dest

            EPQ0(dbg > 6, "orig seq (%lu, %lu), %u\t", l, dest - kc->ts, t);

            *dest |= x & 0xff;
            x >>= 8;

            print_dna(*dest, dbg > 6, '.', 4);

        } else {
            l &= 0xffffffff;             // get length
            ASSERT(l != 0, return -EFAULT);
            uint64_t *s = (uint64_t *)src[0];
            src[1] = (l << BIG_SHFT) | offs;
            dest = kc->ts + (offs >> 2);
            offs += l;
            x = ((unsigned __int128)s[1] << 64) | *s;

            EPQ0(dbg > 6, "orig seqndx (%lu, %lu), %u\t", l, dest - kc->ts, t);
            print_2dna(*dest, x, dbg > 6, 4);
            t <<= 1;
            *dest |= (x << t) & 0xff;
            x >>= 8 - t;

            print_dna(*dest, dbg > 6, '.', 4);

            while (l > 64) {
                *++dest = x & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);

                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 6, '.', 4);

                *++dest = x >>= 8; // zero if t == 0
                print_dna(*dest, dbg > 6, '\n', 4);
                s += 2;
                x = ((unsigned __int128)s[1] << 64) | *s;
                *dest |= (x << t) & 0xff;
                x >>= 8 - t;
                l -= 64;
            }
            print_dna(x, dbg > 6, '.', 4);
            free((uint64_t *)src[0]);
	}
        while (l >= 4) { // first is already written.
            *++dest = x & 0xff;
            print_dna(*dest, dbg > 6, '.', 4);
            x >>= 8;
            l -= 4;
        }
        ++dest;
        ASSERT(dest < kc->ts + ts_end, return -EFAULT, "\n%lu(%lu) > %lu(%lu)?\nkct:0x%x", dest - kc->ts, offs, ts_end, kc->ts_l, kc->ndxkct[ndx]);
        *dest = x; // if l == 0, pre-empty next.
        src[0] = ndx;
        //print_dna(*dest, dbg > 6, '\n', 4);
//EPR0("copied final rest:\t"); print_dna(*dest);
    }
    // why not true, boundaries?
    EPQ(dbg > 4, "offs:%lu kc->ts_l:%lu", offs, kc->ts_l);
    //ASSERT(offs == kc->ts_l, return -EFAULT, "%u, %u", offs, kc->ts_l);
    kc->ts_l = offs;
    // TODO: convert kctndx to nextnt index here.

    return 0;
}

int
fa_read(struct seqb2_t* seq, kct_t* kc)
{
    struct gzfh_t* fhio[3] = { seq->fh, seq->fh + 1, seq->fh + 3};
    int res = -ENOMEM;

    kc->kct = _buf_init_err(kc->kct, 16, goto err);
    kc->id = _buf_init_err(kc->id, 5, goto err);
    kc->s = _buf_init_err(kc->s, 8, goto err);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->ndxkct[i] = -2u;

    /* TODO: load dbSNP and known sites, and mark them. */
    _ACTION(fa_kc(kc, seq->fh + 2), "read and intialized keycounts");
    _ACTION(kct_convert(kc), "converting keycounts");
    _ACTION(save_seqb2(fhio[0], kc), "writing seqb2: %s", fhio[0]->name);
    _ACTION(save_nextnts(fhio[1], kc), "writing next Nts file: %s", fhio[1]->name);
    _ACTION(save_kc(fhio[2], kc), "writing keycounts file: %s", fhio[2]->name);
    _ACTION(reopen(fhio[1], ".nn", ".bd"), "")
    _ACTION(save_boundaries(fhio[1], kc), "writing boundaries: %s", fhio[1]->name);

    res = 0;
err:
    //fclose occurs in main()
    return res;
}

