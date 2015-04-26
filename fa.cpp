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
//#include <pcre.h>
#include "fa.h"


static inline int
print_dna(uint64_t dna, bool dbg = true)
{
    if (dbg == false) return -1;
    for (unsigned t = KEY_WIDTH; t--; dna >>= 2)
        fputc(b6((dna & 3) << 1), stderr);
    fputc('\n', stderr);
    return -1;
}
static inline int
print_dnarc(uint64_t dna, uint64_t rc, bool dbg = true)
{
    if (dbg == false) return -1;
    for (unsigned t = KEY_WIDTH; t--; dna >>= 2)
        fputc(b6((dna & 3) << 1), stderr);
    fputc('|', stderr);
    for (unsigned t = KEY_WIDTH; t--; rc >>= 2)
        fputc(b6((rc & 3) << 1), stderr);
    fputc('\n', stderr);
    return -1;
}

/* Reverse Complement, is ok. */
static inline uint64_t
revcmp(uint64_t dna)
{
    uint64_t m = 0x3333333333333333UL;
    dna = ((dna & m) << 2) | ((dna >> 2) & m);
    m = 0x0f0f0f0f0f0f0f0fUL;
    dna = ((dna & m) << 4) | ((dna >> 4) & m);
    asm ("bswap %0" : "=r" (dna) : "0" (dna));
    dna ^= 0xaaaaaaaaaaaaaaaaUL;
    return dna >> (64 - (KEY_WIDTH << 1));
}

static void
show_list(kct_t* kc, std::list<uint32_t> &bnd)
{
    std::list<uint32_t>::iterator it = bnd.begin();
    unsigned i = 0;
    for (it = bnd.begin(); it != bnd.end(); ++it) {
        Bnd* bd = &kc->bd[*it];
        EPR0("[%u]:\t%u-%u\t%u\t", i++, bd->s, bd->s + bd->l, bd->l);
        print_dnarc(bd->at_dna, bd->dna);
    }
}

static inline int
print_ndx(uint64_t dna, bool dbg = true)
{
    if (dbg == false) return -1;
    uint64_t rc = dna & KEYNT_TRUNC_UPPER;
    dna ^= rc ^ (rc << 1) ^ KEYNT_STRAND;
    rc = revcmp(dna);
    for (unsigned t = KEY_WIDTH; t--; dna >>= 2)
        fputc(b6((dna & 3) << 1), stderr);
    fputc('|', stderr);
    for (unsigned t = KEY_WIDTH; t--; rc >>= 2)
        fputc(b6((rc & 3) << 1), stderr);
    fputc('\n', stderr);
    return -1;
}
enum ensembl_parts {ID, SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META, UNKNOWN_HDR};

enum bnd_type {NEW_REF_ID, N_STRETCH, UQ_REGION, END_REF, DBSNP, KNOWN_SITE};

//ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
// fai does not handle chromosomes with offset.
static Hdr*
new_header(kct_t* kc, void* g, int (*gc) (void*))
{
    int c;
    Hdr* h = new Hdr;
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
                if (++p < ARRAY_SIZE(h->part))
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
                if (++p < ARRAY_SIZE(h->part))
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
    _buf_grow_err(kc->bd, 1ul, 0, return NULL);
    ASSIGN_BD(kc->bd[kc->bd_l], NEW_REF_ID, UNINITIALIZED, 0, UNINITIALIZED, 0, ~0u);

    if (p == NR || p == META) {
        h->p_l = p;
        // ensembl coords are 1-based, we use 0-based.
        kc->bd[kc->bd_l].s = KEY_WIDTH + atoi(kc->id + h->part[START]) - 1;
        h->end_pos = atoi(kc->id + h->part[END]);
        kc->bd[kc->bd_l].l = 0;
    } else {
        h->p_l = UNKNOWN_HDR;
        EPR("\nWARNING: non-ensembl reference.");
    }
    return h;
}

static void
free_kc(kct_t* kc)
{
    for (Kct* k = kc->kct; k != &kc->kct[kc->kct_l]; ++k)
        if (k->seq.m > 3)
            free(k->p.b2);

    std::map<char*, Hdr*, Tid>::iterator it;
    for (it = kc->hdr.begin(); it != kc->hdr.end(); ++it) {
        free(it->second->s);
        delete it->second;
    }
    _buf_free(kc->bd);
    _buf_free(kc->id);
    _buf_free(kc->kct);
    _buf_free(kc->wlkr);
    _buf_free(kc->wbuf);
    _buf_free(kc->kcsndx);
}

/*
 * proh.p_lcess fasta, store 2bit string and count occurances and locations of keys.
 * store per key the count and the last position.
 */
static int
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
        if (c != '>') { // N-stretch
            t = KEY_WIDTH; // use  to count stretch (key is rebuilt later)
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
            ungc(c, g);
            _buf_grow_err(kc->bd, 1ul, 0, return -ENOMEM);
            ASSIGN_BD(kc->bd[kc->bd_l], N_STRETCH, UNINITIALIZED, pos, t, 0, dna);
        } else { // header
            if (h != NULL) {
                // add enddna as length to first element, the chromome boundary.
                //DEBUG_ASSIGN_ENDDNA(kc->bd[*(h->bnd.begin())].at_dna, dna);
                _buf_grow_err(kc->bd, 1ul, 0, return -ENOMEM);
                ASSIGN_BD(kc->bd[kc->bd_l], END_REF, UNINITIALIZED, pos, 0, 0, dna);
                h->bnd.push_back(kc->bd_l++);
            }
            h = new_header(kc, g, gc);
            if (h == NULL) return -EFAULT;
        }
        for (t = 0; t != KEY_WIDTH; ++t) {
            while (isspace(c = gc(g)));
            b ^= b;
            switch(c) {
                case 'C': case 'c': b = 2;
                case 'G': case 'g': b ^= 1;
                case 'U': case 'u': case 'T': case 't': b ^= 2;
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
        _buf_grow_err(kc->bd, 1ul, 0, return -ENOMEM);
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
            if (pos != h->end_pos) {
                EPR("pos != h->end_pos: %u == %u (%u)", pos, h->end_pos, kc->bd[t].l);
                print_dna(dna);
                t = *(h->bnd.begin());
                // correct it.
                kc->bd[t].l = pos;
            }
        } else {
            EPR("=>\tN-stretch at Nt %u", pos);
        }
        //DEBUG_ASSIGN_ENDDNA(kc->bd[kc->bd_l].at_dna, dna);
        ++kc->bd_l;
    }
    if (h == NULL)
        return -EFAULT;
    // add enddna as length to first element, the chromome boundary.
    //DEBUG_ASSIGN_ENDDNA(kc->bd[*(h->bnd.begin())].at_dna, dna);
    _buf_grow_err(kc->bd, 1ul, 0, return -ENOMEM);
    ASSIGN_BD(kc->bd[kc->bd_l], END_REF, UNINITIALIZED, pos, 0, 0, dna);
    h->bnd.push_back(kc->bd_l++);
    return 0;
}

static const unsigned long dbgndx = UNINITIALIZED;

// TODO: instead of excising, move to after (current) last instead.
static int
decr_excise(kct_t const *const kc, unsigned i, const unsigned left)
{
    Walker* wlkr = kc->wlkr;
    uint32_t* wbuf = kc->wbuf;
    while (--i > left) {
        //EPR0("%u/%u ", i, left);
        ASSERT(wbuf[i] != UNINITIALIZED, return -EFAULT, "%u/%u", i, left);
        Kct *x = &_get_kct(kc, wbuf[i]);
        Walker* w = &get_w(wlkr, kc, wbuf[i]);
        // excise out twobits, i.e. shift 2bits above current 2bit
        // down. current bit pos is stored in walker.
        uint8_t* q, *qe;
        ASSERT(w->tmp_count > 0, return -EFAULT);
        --w->tmp_count;
        uint64_t t = w->count;
        if (x->seq.m == 3) {
            q = &x->seq.b2[t >> 2];
            qe = &x->seq.b2[(x->seq.l-1) >> 2];
            --x->seq.l;
        } else {
            q = &x->p.b2[t >> 2];
            qe = &x->p.b2[(x->p.l-1) >> 2];
            // could convert to Kct.seq, but why bother?
            --x->p.l;
        }
        // mask to cover this and past nucleotide.
        t = (1 << (((w->count & 3) + 1) << 1)) - 1;
        *q = ((*q & (t ^ 0xff)) >> 2) | (*q & (t >> 2));
        EPQ(wbuf[i] == dbgndx, "=====> excised <======, became %x", *q | (q != qe ? (q[1] & 3) << 6 : 0));
        wbuf[i] = UNINITIALIZED;
        while (q != qe) {
            *q |= (q[1] & 3) << 6;
            *++q >>= 2;
        }
    }
    //Walker* w = &get_w(wlkr, kc, wbuf[i]);
    //ASSERT(w->tmp_count > 0, return -EFAULT);
    //--w->tmp_count;
    //EPQ(wbuf[i] == dbgndx, "=====> after excised <======");
    wbuf[i] = UNINITIALIZED;

    while (i != 0) {
        if (wbuf[i] != UNINITIALIZED) {
            Walker* w = &get_w(wlkr, kc, wbuf[i]);
            ASSERT(wbuf[i] == UNINITIALIZED, return print_ndx(wbuf[i]), "index %u/%u:%u", i, left, w->tmp_count);
        }
        --i;
    }

    return left;
}

static int dbg = 0;
static unsigned iter = 0;

static inline void merge(Bnd *dest, Bnd *next)
{
    EPQ(dbg > 1, "Extending %u-%u til %u", dest->s, next->s, next->s + next->l);
    dest->l += next->l;
    dest->dna = next->dna;
}

// XXX: inferiority ok? XXX: store deviant bit when unique,
// reverse seq could be required for mapping - maybe 
static int
ext_uq_bnd(kct_t* kc, Hdr* h, Bnd *last)
{
    // N-stretches, later also skips, future: SNVs, splice sites?
    int uqct = 0;
    uint32_t pos = last->s + last->l;
    Bnd *next = &kc->bd[*h->bdit];
    const char* dbgtid = "GL000210.1 ";//"GL000207.1";//MT"; //GL000239.1";
    const char* hdr = kc->id + h->part[0];
    dbg = (strncmp(hdr, dbgtid, strlen(dbgtid)) == 0) * 2;

    uint64_t t, dna = last->dna; // first seq after skip
    uint64_t rc = revcmp(dna);
    uint64_t ndx = _getxtdndx(kc, ndx, dna, rc);
    Walker* w = &get_w(kc->wlkr, kc, ndx);
    uint32_t infior = max(last->i + 1, w->infior); //XXX
    if (dbg > 2)
        show_list(kc, h->bnd);

    _buf_grow0(kc->bd, 1ul);
    Bnd *inter = &kc->bd[kc->bd_l];
    ASSIGN_BD(*inter, last->t, UNINITIALIZED, pos, 0, infior, dna);

    // boundary is considered as a first unique,
    // if we find a 2nd, we really want to extend this region
    const int ext = kc->ext;
    int left = ext;
    EPQ0(dbg > 0, "----[\t%s%s:(%u-%u)..(%u-%u)\t(%lu)@%u\t]----\tdna:",
            strlen(hdr) < 8 ? "\t":"",hdr, last->s, last->s+last->l, next->s, next->s+next->l, last->t, pos); print_dna(dna, dbg > 0);

    for (;pos < next->s;++pos) { // until next event
        EPQ(ndx == dbgndx, "observed dbgndx at %u", pos);
        ASSERT(_getxtdndx0(kc, ndx) < kc->kct_l, return -print_dna(dna), "at %u, 0x%lx", pos, ndx);
        EPQ0(dbg > 2, "%u:\t", pos);
        print_dna(dna, dbg > 2);

        w = &get_w(kc->wlkr, kc, ndx);
        if (left) { // keep filling the buffer
            if (left == ext) {
                //EPR("position after 2nd unique, or each extended.");
                // what we'll jump to, in subsequent iterations.
                // This is inserted (once) when region is complete (left becomes 0) -
                // otherwise it's overwritten successively.
                inter->dna = dna;
                inter->l = pos - inter->s;
            }
            if (--left) {
                ASSERT(kc->wbuf[left] == UNINITIALIZED, return -EFAULT, "[%u]", pos & print_ndx(ndx));
                kc->wbuf[left] = ndx;
                if (infior > w->infior)
                    w->infior = infior;
            } else { // actually buf just became full, without extension(s not already handled)..
                EPQ(inter->s > 0xfffffff, "%u", -inter->s);
                infior = 0;

                if (inter->l != 0) {
                    // we did insert a range
                    EPQ0(dbg > 1, "\n[%u-%u]\t", inter->s, inter->s + inter->l - 1);print_dnarc(inter->at_dna, inter->dna, dbg > 1);
                    //show_list(kc, h->bnd);
                    if ((last->s + last->l) < inter->s) {
                        _buf_grow0(kc->bd, 2ul);
                        h->bnd.insert(h->bdit, kc->bd_l); // even after increment one should be available for inter.
                        last = kc->bd + kc->bd_l++;
                        inter = last + 1;
                    } else { // join - may be undesirable for certain boundary types in the future
                        merge(last, inter);
                    }
                    next = &kc->bd[*h->bdit];
                    inter->l = 0;
                }
                // add tmp_count for each - we cannot extend it in this iteration
                for (left = ext; --left;) {
                    get_w(kc->wlkr, kc, kc->wbuf[left]).count++;
                    --get_w(kc->wlkr, kc, kc->wbuf[left]).tmp_count;
                    kc->wbuf[left] = UNINITIALIZED;
                }
            }
        }
        t = w->count + w->tmp_count;
        Kct *y = &_get_kct(kc, ndx);
        uint8_t b2;
        if (y->seq.m == 3 || y->seq.m == 2) { //EPR("Kct in seq format");
            b2 = (y->seq.b2[t >> 2] >> ((t & 3) << 1)) & 3;
            t = y->seq.l;
        } else {
            b2 = (y->p.b2[t >> 2] >> ((t & 3) << 1)) & 3; // get next 2bit
            t = y->p.l; // can also be 1 (unique), unless we decide to convert back upon decrement
        }
        if (h->s) {
            uint8_t sb2 = (h->s[pos>>2] >> ((pos & 3) << 1)) & 3;
            ASSERT(b2 == sb2, return print_dnarc(dna, rc), "[%u]: expected %c, got %c for ndx 0x%lx", pos,b6(sb2<<1), b6(b2<<1), ndx);
        }
        dna = _seq_next(b2, dna, rc);
        ndx = _getxtdndx(kc, ndx, dna, rc);
        if (t > 1ul) {
            if (left == 0) ++w->count;
            else ++w->tmp_count;
            continue;
        }
        EPQ0(dbg > 1, "%u,", pos);
        // found unique, store orientation of 2nd bit of central Nt.
        // if first bit is set, the key is in the wrong orientation.
        // in any case - before any revcmping, the central bit needs to be
        // put back in to ge the sequence at this site.
        y->seq.m ^= (y->seq.m & 1) ^ !(dna & KEYNT_STRAND);
        ++infior;
        if (left == 0) { // this is a first unique.
            ASSIGN_BD(*inter, UQ_REGION, UNINITIALIZED, pos, 0, max(infior, w->infior + 1u), dna);
            left = ext;
            infior = inter->i + 1;
            continue;
        }
        // insert a new range.
        ++uqct;
        left = decr_excise(kc, ext, left);
        if (left < 0) return left;
        left = ext;
    }
    ASSERT(pos == next->s, return -EFAULT);
    EPQ(dbg > 2, "---");
    //ASSERT(dna == next->at_dna, return -print_dnarc(dna, next->at_dna), "[%u] => enddna mismatch", pos);
    if (left) {
        ASSERT((inter->s + inter->l + ext) >= next->s, return -EFAULT);
        EPQ0(dbg > 1, "\nExtending %u-%u til %u\t(next:%u-%u)\t", last->s, last->s+last->l, inter->s+inter->l, next->s, next->s+next->l); print_dnarc(last->dna, inter->at_dna, dbg > 1);
        next->l += next->s - inter->s;
        next->s = inter->s;
        next->at_dna = inter->at_dna;
        //XXX: why decrement left?
        left = decr_excise(kc, ext, --left);
        if (left < 0)
            return left;
        // if last and next have now become adjoining, merge the two
        // XXX why not (last->s + last->l + ext > next->s) ? (assertion error later, but why?)
        if ((last->s + last->l >= next->s)) {
            if (next->t != UQ_REGION)
                last->t = next->t;
            EPQ(dbg > 1,"merging %u-%u", last->s, next->s + next->l);
            merge(last, next);
            h->bdit = h->bnd.erase(h->bdit);
            --h->bdit;
        }

    } else {
        EPQ(dbg > 1, "Not joining\n");
    }
//                pos += kc->bd[*(h->bdit)].l; // add skipped Nts
    return uqct;
}

static int
ext_uq_hdr(kct_t* kc, Hdr* h)
{
    int uqct = 0;
    h->bdit = h->bnd.begin();

    for (Bnd *last = &kc->bd[*h->bdit]; ++h->bdit != h->bnd.end(); last = &kc->bd[*h->bdit]) {
        int ret = ext_uq_bnd(kc, h, last);
        if (ret < 0) return ret;
        uqct += ret;
    }
    for (unsigned t = 0; t != kc->kct_l; ++t) { //XXX
        Walker* w = kc->wlkr + t;
        //XXX
        ASSERT(w->tmp_count == 0u, return -EFAULT, "%u/%u", t, kc->kct_l);
    }
    if (dbg > 1)
        show_list(kc, h->bnd);

    return uqct;
}

static int
ext_uq_iter(kct_t* kc)
{
    int uqct = 0;
    for (std::list<Hdr*>::iterator h = kc->h.begin(); h != kc->h.end(); ++h) {
        // over ref headers
        int ret = ext_uq_hdr(kc, *h);
        if (ret < 0) return ret;
        uqct += ret;
    }
    EPR("extended %u unique ranges in iteration %u", uqct, iter++);
    Walker * w = kc->wlkr;
    for (unsigned t = 0; t != kc->kct_l; ++t, ++w) {
        ASSERT(w->tmp_count == 0u, return -EFAULT, "%u/%u", t, kc->kct_l);
        w->count = 0u;
    }
    return uqct;
}

int
fn_convert(struct gzfh_t* fhio, const char* search, const char* replace)
{
    char* f = strstr(fhio->name, search);
    if (f == NULL) return 0;
    strncpy(f, replace, strlen(replace) + 1);
    return 1;
}

int write1(struct gzfh_t* fhout, kct_t* kc)
{
    int ret = -EFAULT;
    if (fhout->write(fhout, (const char*)&kc->bd_l, sizeof(uint32_t)) < 0)
        return ret;
    EPR("%u", kc->bd_l);
    if (fhout->write(fhout, (const char*)&kc->id_l, sizeof(uint32_t)) < 0)
        return ret;
    if (fhout->write(fhout, (const char*)kc->id, kc->id_l * sizeof(char)) < 0)
        return ret;
    if (fhout->write(fhout, (const char*)&kc->kct_l, sizeof(uint32_t)) < 0)
        return ret;
    uint32_t len = kc->h.size();
    if (fhout->write(fhout, (const char*)&len, sizeof(uint32_t)) < 0)
        return ret;
    for(std::list<Hdr*>::iterator h = kc->h.begin(); h != kc->h.end(); ++h)
    {
        EPQ(dbg > 1, "%s..", kc->id + (*h)->part[0]);
        if (fhout->write(fhout, (const char*)&(*h)->end_pos, sizeof(uint32_t)) < 0)
            return ret;
        uint32_t len = (*h)->bnd.size();
        if (fhout->write(fhout, (const char*)&len, sizeof(uint32_t)) < 0)
            return ret;
        for(std::list<uint32_t>::iterator b = (*h)->bnd.begin(); b != (*h)->bnd.end(); ++b)
        {
            if (fhout->write(fhout, (const char*)&*b, sizeof(uint32_t)) < 0)
                return ret;
        }
        // sequence is not stored.
        if (fhout->write(fhout, (const char*)(*h)->part, 10 * sizeof(uint16_t)) < 0)
            return ret;
        if (fhout->write(fhout, (const char*)&(*h)->p_l, sizeof(uint8_t)) < 0)
            return ret;
        len = (*h)->bnd.size();
        // p_l is not needed anymore
    }
    if (fhout->write(fhout, (const char*)&kc->bd_m, sizeof(uint8_t)) < 0)
        return ret;
    // FIXME: full flush silences valgrind errors on close, about uninit
    // values - either when writing buf or kc->kct.
    ASSERT(gzflush(fhout->io, Z_FULL_FLUSH) == Z_OK, return -EFAULT);
    EPQ(dbg > 0, "after flush1");
    if (fhout->write(fhout, (const char*)kc->kct, kc->kct_l * sizeof(Kct)) < 0)
        return ret;


    len = kc->kct_l * 2;
    uint32_t* buf = (uint32_t*)malloc(sizeof(uint32_t) * len);
    if (buf == NULL) return -ENOMEM;
    uint32_t i = 0;
    // FIXME: may need 64 bit index.
    for (uint32_t ndx = 0u; i != len; ++ndx) {
        if (kc->kcsndx[ndx] >= kc->kct_l) {
            kc->kcsndx[ndx] = kc->kct_l; // make high bit range available
        } else {
            buf[i++] = ndx;
            buf[i++] = kc->kcsndx[ndx];
        }
    }
    if (fhout->write(fhout, (const char*)buf, sizeof(uint32_t) * len) < 0)
        goto err;
    // kc->wlkr and kc->wbuf aren't filled yet.
    if (fhout->write(fhout, (const char*)kc->bd, kc->bd_l * sizeof(Bnd)) < 0)
        return ret;
    ret = 0;
err:
    free(buf);
    return ret;
}

int restore1(struct gzfh_t* fhin, kct_t* kc)
{
    int ret = -EFAULT;
    if (fhin->read(fhin, (char*)&kc->bd_l, sizeof(uint32_t)) < 0)
        return ret;
    EPR("%u", kc->bd_l);
    if (fhin->read(fhin, (char*)&kc->id_l, sizeof(uint32_t)) < 0)
        return ret;
    uint32_t len = kc->id_l * sizeof(char);
    kc->id = (char*)malloc(len);
    if (fhin->read(fhin, (char*)kc->id, len) < 0)
        return ret;
    if (fhin->read(fhin, (char*)&kc->kct_l, sizeof(uint32_t)) < 0)
        return ret;
    uint32_t blen, val;
    if (fhin->read(fhin, (char*)&len, sizeof(uint32_t)) < 0)
        return ret;
    for (uint32_t i=0; i != len; ++i) {
        Hdr *h = new Hdr;
        h->s = NULL; h->s_l = 0u; // XXX: may be removed in future.
        if (fhin->read(fhin, (char*)&h->end_pos, sizeof(uint32_t)) < 0)
            return ret;
        if (fhin->read(fhin, (char*)&blen, sizeof(uint32_t)) < 0)
            return ret;
        for (uint32_t j=0; j != blen; ++j) {
            if (fhin->read(fhin, (char*)&val, sizeof(uint32_t)) < 0)
                return ret;
            h->bnd.push_back(val);
        }
        if (fhin->read(fhin, (char*)h->part, 10 * sizeof(uint16_t)) < 0)
            return ret;
        if (fhin->read(fhin, (char*)&h->p_l, sizeof(uint8_t)) < 0)
            return ret;
        kc->h.push_back(h);
        kc->hdr.insert(std::pair<char*, Hdr*>(kc->id + h->part[ID], kc->h.back()));
    }
    if (fhin->read(fhin, (char*)&kc->bd_m, sizeof(uint8_t)) < 0)
        return ret;
    len = kc->kct_l * sizeof(Kct);
    kc->kct = (Kct*)malloc(len);
    if (fhin->read(fhin, (char*)kc->kct, len) < 0)
        return ret;
    memset(kc->kcsndx, kc->kct_l, KEYNT_BUFSZ * sizeof(kc->kcsndx[0]));
    for (uint32_t i=0; i != kc->kct_l; ++i) {
        if (fhin->read(fhin, (char*)&val, sizeof(uint32_t)) < 0)
            return ret;
        if (fhin->read(fhin, (char*)&len, sizeof(uint32_t)) < 0)
            return ret;
        kc->kcsndx[val] = len;
    }
    len = kc->bd_l * sizeof(Bnd);
    kc->bd = (Bnd*)malloc(len);
    if (fhin->read(fhin, (char*)kc->bd, len) < 0)
        return ret;
    return 0;
}

int fa_index(struct seqb2_t* seq, uint32_t blocksize)
{
    struct gzfh_t* fhio = seq->fh;
    struct gzfh_t* fhin = seq->fh + 2; // init with ref
    size_t t;
    int res, ret = -1;
    void* g;
    int (*gc) (void*);
    int (*ungc) (int, void*);
    int is_gzfile = fhin->io != NULL;
    char file[256];
    kct_t kc;
    kc.ext = seq->readlength - KEY_WIDTH;
    kc.kcsndx = _buf_init_arr(kc.kcsndx, KEYNT_BUFSZ_SHFT);

    const char* ndxact[4] = {".fa", "_x1.gz"};

    if (fhio->name != NULL && ((res = strlen(fhio->name)) > 255)) {
        strncpy(file, fhio->name, res);
    } else {
        res = strlen(fhin->name);
        ASSERT(strncmp(fhin->name, "stdin", res) != 0, return -1);
        strcpy(file, fhin->name);
    }
    fhio->name = &file[0];
    if (!fn_convert(fhio, ndxact[0], ndxact[1]))
        return -1; // skip if this file was not provided on the commandline

    fprintf(stderr, "== %s(%lu)\n", fhio->name, fhio - seq->fh);
    res = set_io_fh(fhio, blocksize, 0);
    if (res < 0) { fprintf(stderr, "== set_io_fh failed:%d", res); goto out; }

    if (res == 1) {
        kc.kct = _buf_init(kc.kct, 16);
        kc.bd = _buf_init(kc.bd, 1);

        // FIXME: assigning 1 Mb for reference id, to keep realloc from happening.
        // (realloc would invalidate pointers inserted in kc.hdr)
        kc.id = _buf_init(kc.id, 20);

        memset(kc.kcsndx, UNINITIALIZED, KEYNT_BUFSZ * sizeof(kc.kcsndx[0]));

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
        res = fa_kc(&kc, g, gc, ungc);
        EPR("left fa_kc");fflush(NULL);
        ASSERT(res >= 0, goto out);

        // write to disk, TODO: load from disk.
        res = write1(fhio, &kc);
        //ASSERT(res >= 0, goto out);
        //res = fhio->close(fhio->io);
    } else {
        res = restore1(fhio, &kc);
        ASSERT(res >= 0, goto out);
        res = fhio->close(fhio->io);
    }
    ASSERT(res >= 0, goto out);

    t = kc.kct_l * sizeof(Walker);
    kc.wlkr = (Walker*)malloc(t);
    memset(kc.wlkr, 0ul, t);

    t = kc.ext * sizeof(uint32_t);
    kc.wbuf = (uint32_t*)malloc(t);
    memset(kc.wbuf, UNINITIALIZED, t);

    do { // until no no more new uniques
        res = ext_uq_iter(&kc);
    } while (res > 0);
    ASSERT(res >= 0, goto out);

    ret = res != 0;

out:
    free_kc(&kc);
    return ret;

}


