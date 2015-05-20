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

    if (p == NR || p == META) {
        h->p_l = p;
        // ensembl coords are 1-based, we use 0-based.
        kc->bd[kc->bd_l].l = KEY_WIDTH + atoi(kc->id + h->part[START]) - 1
            -kc->bd[kc->bd_l].s;
        h->end_pos = atoi(kc->id + h->part[END]);
    } else { //TODO: hash lookup from fai
        kc->bd[kc->bd_l].l = -kc->bd[kc->bd_l].s;
        h->p_l = UNKNOWN_HDR;
        EPR("\nWARNING: non-ensembl reference.");
    }
    return h;
}
/*
 * process fasta, store 2bit string and count occurances and locations of keys.
 * store per key the count and the last position.
 */
static int
fa_kc(kct_t* kc, void* g, int (*gc) (void*), int (*ungc) (int, void*))
{
    uint64_t dna = 0, rc = 0, ndx, b = 0, t = 0ul;
    uint32_t pos = 0u;
    int c;
    Hdr* h = NULL;
    kc->s_l = 0ul;
    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID

    // TODO: realpos based dbsnp or known site boundary insertion.
    while (c >= 0) {
        kc->bd[kc->bd_l].s = pos;
        if (c != '>') { // N-stretch
            do {
                if (c != '\n') {
                    if (c != 'N') {
                        if (B6(b, c) || c == '>') break;
                        EPR("Strange nucleotide:%c, treated as N\n", c);
                    }
                    _addtoseq(kc->s, 0);
                    ++t;
                }
                c = gc(g);
            } while(c != -1);
            // TODO: handle single or a few N's differently.

            EPR("=>\tN-stretch at Nt %u (%lu)", pos, t - KEY_WIDTH);
            if (c == '>' || c == -1) { // skip N's at end.
                EPR("processed %u Nts for %s", pos, kc->id + h->part[0]);
                // pos += t - KEY_WIDTH;
                // non-fatal for Y
                EPQ(pos + t - KEY_WIDTH != h->end_pos, "pos + t - KEY_WIDTH != h->end_pos: %u (+'ed:%lu) == %u",
                        pos, t, h->end_pos);
                t = 0;
                continue;
            }
            kc->bd[kc->bd_l].l = t;
        } else {
            // Append the index of the next boundary to the last header.
            // The at_dna is of the last chromo. XXX if at_dna is removed
            // altogether we can stop at h->bnd.end() instead.
            if (h != NULL)
                h->bnd.push_back(kc->bd_l);

            h = new_header(kc, g, gc);
            h->s_s = kc->s_l; 
            if (h == NULL) return -EFAULT;
            EPQ(dbg > 2, "header %s", kc->id + h->part[0]);
            c = gc(g);
            if (!isb6(b6(c))) {
                t = -pos + KEY_WIDTH;
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
default: // include 'N's and such to make sure the key is completed.
                    _addtoseq(kc->s, b);
                    dna = _seq_next(b, dna, rc);
                    continue;
            }
            break;
        }
        kc->bd[kc->bd_l].i = 0;
        pos = kc->bd[kc->bd_l].s + kc->bd[kc->bd_l].l;
        if_ever (t != KEY_WIDTH) { // c == -1, 'N' or '>'
            // we can get here if an N-stretch follows a new contig
            if (t != 0) {
                // or if two N-s(tretches) lie close to one anorther
                EPQ(dbg > 1, "Key incomplete before next boundary - merging");
            }
            EPQ(dbg > 1, "jump at %u", pos);
            continue; // join regions.
        }
        // ndx (key) is 2nd cNt bit excised, rc or dna
        kc->bd[kc->bd_l].dna = dna;
        h->bnd.push_back(kc->bd_l);

        while (1) {
            if (isspace(c = gc(g))) continue;
            uint64_t *ct;
            _get_ndx(ndx, ndx, dna, rc);
            if (kc->kctndx[ndx] != ~0ul) {
                ASSERT(kc->kctndx[ndx] < kc->kct_l, return -EFAULT, "0x%lx, 0x%lx", kc->kctndx[ndx], ndx);
                ct = kc->kct + kc->kctndx[ndx];
                // TODO: using a length - based conversion, we could cram in 30th bit.
                if (*ct & FLAG_B2CT) {
                    EPQ(dbg > 3, "sequence & length");
                    t = *ct >> B2LEN_OFFS_SHFT;

                    if (t < 29ul) {
                        *ct += B2LEN_OFFS;
                    } else {
                        EPQ(dbg > 3, "conversion to index:%c", c);
                        t = *ct;                                  // temp store 2bit dna
                        *ct = ndx = kc->kce_l;                    // conversion: set to index.
                        ct = (uint64_t*)malloc(sizeof(uint64_t)); // make room for 2bit dna
                        _buf_grow0(kc->kce, 1ul);
                        kc->kce[kc->kce_l++] = { .m = 0, .l = 30, .b2 = ct }; // 29th will be added
                        *ct = t & (FLAG_B2CT - 1ul);     // after current 2bit dna, in ndx.
                        t = 29;
                    }
                } else {
                    EPQ(dbg > 3, "in index format:%c", c);
                    ndx = *ct & INDEX_MASK;

                    t = kc->kce[ndx].l++;                // *ct is index to an extended keycount
                    if (t & 0x1f) {
                        ct = kc->kce[ndx].b2 + (t >> 5); // point to the uint32_t to be updated
                    } else {
                        if (t == (0x20u << kc->kce[ndx].m)) {
                            kc->kce[ndx].b2 = (uint64_t *)realloc(kc->kce[ndx].b2,
                                    8 << ++kc->kce[ndx].m);
                            ASSERT(kc->kce[ndx].b2 != NULL, return -ENOMEM);
                        }
                        ct = kc->kce[ndx].b2 + (t >> 5);
                        *ct = 0ul;
                    }
                    t &= 0x1f;
                }
            } else {
                EPQ(dbg > 3, "new key 0x%lx at %u", ndx, pos);
                _buf_grow(kc->kct, 1ul, 0);
                kc->kctndx[ndx] = kc->kct_l;
                ct = kc->kct + kc->kct_l++;
                *ct = B2LEN_OFFS | FLAG_B2CT;
                t = 0;
            }
            b ^= b;
            switch (c) {
case 'C': case 'c': b = 2;
case 'G': case 'g': b ^= 1;
case 'T': case 't':
case 'U': case 'u': b ^= 2;
case 'A': case 'a':
                // append next 2bit to last ndx. Early Nts are in low bits.
                *ct |= b << (t << 1);
                EPQ0(dbg > 3, "%s\t%u\t", kc->id + h->part[0], pos);
                print_dna(*ct, dbg > 3, '\n', (ct >= kc->kct) && (ct < (kc->kct + kc->kct_l)) ? 29 : 32);
                ++pos;
                _addtoseq(kc->s, b);
                dna = _seq_next(b, dna, rc);
                continue;
            }
            // happens at a boundary: no next Nts.
            if ((ct >= kc->kct) && (ct < (kc->kct + kc->kct_l)))
                *ct -= B2LEN_OFFS;
            else // If not within ->kct range then this is an extended keycount:
                --kc->kce[ndx].l;
            break;
        }

        kc->bd[kc->bd_l++].at_dna = dna;
        _buf_grow0(kc->bd, 2ul);
        if (c == '>' || c == -1) {
            EPR("processed %u Nts for %s", pos, kc->id + h->part[0]);
            // non-fatal for Y
            EPQ(pos != h->end_pos, "pos - KEY_WIDTH != h->end_pos: %u == %u",
                    pos, h->end_pos);
        }
        t = KEY_WIDTH;
    }
    if (h == NULL)
        return -EFAULT;
    kc->bd[0].l -= pos;
    kc->bd[0].s = pos;
    h->bnd.push_back(0);
    return 0;
}

// concatenate 2bits, byte aligned, for each key in buffer kc->ts,
// store lengths in kc->kct;
static int
kct_convert(kct_t* kc)
{
    kc->ts_l = kc->s_l; // about as many next NTs as Nts.

    const uint64_t ts_end = (kc->ts_l >> 2) + !!(kc->ts_l & 3);
    uint8_t *dest = (uint8_t*)malloc(ts_end);
    ASSERT(dest != NULL, return -ENOMEM);
    kc->ts = dest;
    *dest = '\0';
    uint64_t offs = 0ul;

    for (uint64_t *src = kc->kct; src != &kc->kct[kc->kct_l]; ++src) {
	uint64_t l, x = *src;
        unsigned t = offs & 3;
        ASSERT(dest < kc->ts + ts_end, return -EFAULT);
        if (x & FLAG_B2CT) { // then x contains count and next Nts below this bit

            l = x >> B2LEN_OFFS_SHFT;    // isolate twobit count
            *src = offs + l;             // kc->kct conversion: store end position
                                         // of its next Nts from now on.
            dest = kc->ts + (offs >> 2); // destination for copy.
            x &= FLAG_B2CT - 1ul;        // remove len => what is left is 2bit.
            x <<= t << 1;                // prealign 2bits for appending to dest
            offs += l;                   // update offs for next round.

            dbg = offs == dbgoffs ?  dbg | 4 : dbg & ~4;
            EPQ0(dbg > 3, "orig seq (%lu, %lu), %u\t", l, dest - kc->ts, t);
            print_2dna(*dest, x, dbg > 3, 4);

            *dest |= x & 0xff;
            x >>= 8;

            print_dna(*dest, dbg > 3, '.', 4);

        } else {
            kct_ext* ke = kc->kce + x;   // x was an index to a extended keycount.
            l = ke->l;                   // get length
            *src = offs + l;             // kc->kct conversion
            dest = kc->ts + (offs >> 2);
            uint64_t *s = ke->b2;
            x = *s;

            dbg = (offs + l) == dbgoffs ?  dbg | 4 : dbg & ~4;
            EPQ0(dbg > 3, "orig seqndx (%lu, %lu), %u\t", l, dest - kc->ts, t);
            print_2dna(*dest, x, dbg > 3, 4);
            t <<= 1;
            *dest |= (x << t) & 0xff;
            x >>= 8 - t;
            offs += l;                  // update offs for next.

            print_dna(*dest, dbg > 3, '.', 4);

            while (l > 32) {
                *++dest = x & 0xff;
                print_dna(*dest, dbg > 3, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 3, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 3, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 3, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 3, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 3, '.', 4);
                *++dest = (x >>= 8) & 0xff;
                print_dna(*dest, dbg > 3, '.', 4);
                *++dest = x >>= 8; // zero if t == 0
                print_dna(*dest, dbg > 3, '\n', 4);
                x = *++s;
                *dest |= (x << t) & 0xff;
                x >>= 8 - t;
                l -= 32;
            }
            print_dna(x, dbg > 3, '.', 4);
            free(ke->b2);
	}
        while (l >= 4) { // first is already written.
            *++dest = x & 0xff;
            print_dna(*dest, dbg > 3, '.', 4);
            x >>= 8;
            l -= 4;
        }
        *++dest = x; // if l == 0, pre-empty next.
        print_dna(*dest, dbg > 3, '\n', 4);
//EPR0("copied final rest:\t"); print_dna(*dest);
    }
    // why not true, boundaries?
    //ASSERT(offs == kc->ts_l, return -EFAULT, "%u, %u", offs, kc->ts_l);
    kc->ts_l = offs;
    // TODO: convert kctndx to nextnt index here.

    return 0;
}

int
fa_read(struct seqb2_t* seq, kct_t* kc)
{
    int res = -ENOMEM;
    void* g;
    int (*gc) (void*);
    int (*ungc) (int, void*);
    struct gzfh_t* fhin = seq->fh + 2; // init with ref
    struct gzfh_t* fhio[3] = { seq->fh, seq->fh + 1, seq->fh + 3};
    int is_gzfile = fhin->io != NULL;

    kc->kct = _buf_init_err(kc->kct, 16, goto err);
    kc->kce = _buf_init_err(kc->kce, 8, goto err);
    kc->bd = _buf_init_err(kc->bd, 1, goto err);
    kc->id = _buf_init_err(kc->id, 5, goto err);
    kc->s = _buf_init_err(kc->s, 8, goto err);

    for (uint64_t i=0ul; i != KEYNT_BUFSZ; ++i)
        kc->kctndx[i] = ~0ul;

    if (is_gzfile) {
        g = fhin->io;
        gc = (int (*)(void*))&gzgetc;
        ungc = (int (*)(int, void*))&gzungetc;
    } else {
        g = fhin->fp;
        gc = (int (*)(void*))&fgetc;
        ungc = (int (*)(int, void*))&ungetc;
    }
    /* TODO: load dbSNP and known sites, and mark them. */
    _ACTION(fa_kc(kc, g, gc, ungc), "read and intialized keycounts");
    _ACTION(kct_convert(kc), "converting keycounts");
    _ACTION(save_seqb2(fhio[0], kc), "writing seqb2: %s", fhio[0]->name);
    _ACTION(save_nextnts(fhio[1], kc), "writing next Nts file: %s", fhio[1]->name);
    _ACTION(save_kc(fhio[2], kc), "writing keycounts file: %s", fhio[2]->name);
    _ACTION(reopen(fhio[1], ".nn", ".bd"), "")
    _ACTION(save_boundaries(fhio[1], kc), "writing boundaries: %s", fhio[1]->name);

    res = 0;
err:
    //fclose occurs in main()
    _buf_free(kc->kce);
    return res;
}

