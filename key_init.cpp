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

#define ENS_HDR_PARTCT 10
//ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
// fai does not handle chromosomes with offset.
static Hdr*
new_header(kct_t* kc, void* g, int (*gc) (void*))
{
    int c;
    Hdr* h = new Hdr;
    h->part = (uint16_t*)malloc(ENS_HDR_PARTCT * sizeof(uint16_t));
    ASSERT(h->part != NULL, return NULL);

    h->s = _buf_init_err(h->s, 8, return NULL);
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
                    char* q = &kc->id[h->part[p]];
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
                    char* q = &kc->id[h->part[p]];
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

    kc->hdr.insert(std::pair<char*, Hdr*>(kc->id + h->part[ID], h));
    if (p == NR || p == META) {
        h->p_l = p;
        // ensembl coords are 1-based, we use 0-based.
        kc->bd[kc->bd_l].l = KEY_WIDTH + atoi(kc->id + h->part[START]) - 1
            -kc->bd[kc->bd_l].s;
        h->end_pos = atoi(kc->id + h->part[END]);
    } else {
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
int
fa_kc(kct_t* kc, void* g, int (*gc) (void*), int (*ungc) (int, void*))
{
    uint64_t dna = 0, rc = 0, ndx, b = 0;
    uint32_t t, pos = 0u;
    int c;
    Kct *ct = NULL;
    uint8_t *q = NULL;
    Hdr* h = NULL;

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID

    // TODO: realpos based dbsnp or known site boundary insertion.
    while (c >= 0) {
        kc->bd[kc->bd_l] = {.at_dna = dna, .dna = 0, .s = pos, .l = 0, .i = 0};
        if (c != '>') { // N-stretch
            t = KEY_WIDTH; // stretch length + shift required for new key
            do {
                if (c != '\n') {
                    if (c != 'N') {
                        if (B6(b, c) || c == '>') break;
                        EPR("Ignoring strange nucleotide:%c\n", c);
                    }
                    _addtoseq(h->s, 0);
                    ++t;
                }
                c = gc(g);
            } while(c != -1);
            if (c == -1) break;
            kc->bd[kc->bd_l].l = t;
            ungc(c, g);
        } else { // header
            // Append the index of the next boundary to the last header.
            // The at_dna is of the last chromo. XXX if at_dna is removed
            // altogether we can stop at h->bnd.end() instead.
            if (h != NULL)
                h->bnd.push_back(kc->bd_l);

            h = new_header(kc, g, gc);
            if (h == NULL) return -EFAULT;
        }
        for (t = 0; t != KEY_WIDTH; ++t) {
            while (isspace(c = gc(g)));
            b ^= b;
            switch(c) {
case 'C': case 'c': b = 2;
case 'G': case 'g': b ^= 1;
case 'U': case 'u':
case 'T': case 't': b ^= 2;
case 'A': case 'a':
                    _addtoseq(h->s, b);
                    dna = _seq_next(b, dna, rc);
                    continue;
            }
            EPR("key incomplete before next boundary - "
                    "may cause enddna assertion failure");
            break;
        }
        if_ever (t != KEY_WIDTH)// c == -1, 'N' or '>'
            continue;
        pos = kc->bd[kc->bd_l].s + kc->bd[kc->bd_l].l;
        // ndx (key) is 2nd cNt bit excised, rc or dna
        kc->bd[kc->bd_l].dna = dna;
        h->bnd.push_back(kc->bd_l);

        do {
            ndx = _get_ndx(ndx, dna, rc);
            if (kc->kcsndx[ndx] != UNINITIALIZED) {
                ct = kc->kct + kc->kcsndx[ndx];
                if (ct->seq.m == 3) { //EPR("Kct has seq format");
                    t = ct->seq.l++;
                    q = &ct->seq.b2[t >> 2];
                    if ((t & 3) == 0) {
                        b = ARRAY_SIZE(ct->seq.b2);
                        if (t == (b << 2)) { // convert to p (pointer) format
                            q = (uint8_t *)malloc(1u << ++ct->seq.m);
                            memcpy(q, ct->seq.b2, b); 
                            ASSERT(ct->seq.m == ct->p.m, return -EINVAL,
                                    "struct not aligned as expected\n");
                            ASSERT(ct->seq.l == (ct->p.l & 0xff), return -EINVAL,
                                    "endian test failed\n"); // redundant?
                            ct->p.l &= 0xff;
                            ct->p.b2 = q;
                            q += b;
                        }
                        *q = '\0';
                    }
                } else { //EPR("Kct in p format");
                    t = ct->p.l++;
                    q = &ct->p.b2[t >> 2];
                    if ((t & 3) == 0) {
                        if (t == (4u << ct->p.m)) {
                            q = (uint8_t *)realloc(ct->p.b2, 1 << ++ct->p.m);
                            if_ever (q == NULL) return -ENOMEM;
                            ct->p.b2 = q;
                            q += t >> 2;
                        }
                        *q = '\0';
                    }
                }
            } else { //EPR("new key 0x%lx at %u", ndx, pos);
                _buf_grow(kc->kct, 1ul, 0);
                kc->kcsndx[ndx] = kc->kct_l;
                ct = kc->kct + kc->kct_l++;
                ct->seq = {.m = 3, .l = 1, 0}; // init remaining to zero
                q = ct->seq.b2;
                *q = '\0';
                t = 0;
            }
            // next iteration ...
            while (isspace(c = gc(g)));
            b ^= b;
            switch(c) {
                case 'C': case 'c': b = 2;
                case 'G': case 'g': b ^= 1;
                case 'T': case 'U': case 't': case 'u': b ^= 2;
                    *q |= b << ((t & 3) << 1); // append next 2bit to last ndx
                case 'A': case 'a':
                    _addtoseq(h->s, b);
                    dna = _seq_next(b, dna, rc);
                    ++pos;
                    continue;
            }
            if (ct->seq.m == 3) --ct->seq.l;
            else --ct->p.l;
            break;
        } while (c >= 0);

        if (c == '>' || c == -1) {
            EPR("processed %u Nts for %s", pos, kc->id + h->part[0]);
            ASSERT (pos == h->end_pos, return -EFAULT, "pos != h->end_pos: %u == %u (%u)", pos, h->end_pos, kc->bd[t].l);
        } else {
            EPR("=>\tN-stretch at Nt %u", pos);
        }
        kc->bd[kc->bd_l++].at_dna = dna;
        _buf_grow0(kc->bd, 2ul);
    }
    if (h == NULL)
        return -EFAULT;
    kc->bd[0].l -= pos;
    kc->bd[0].s = pos;
    h->bnd.push_back(0);
    return 0;
}


