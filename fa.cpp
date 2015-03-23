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
#include <pcre.h>
#include "fa.h"

static inline int print_dna(uint64_t dna, bool dbg = true)
{
    if (dbg == false) return -1;
    for (unsigned t = KEY_WIDTH; t--; dna >>= 2)
        fputc(b6((dna & 3) << 1), stderr);
    fputc('\n', stderr);
    return -1;
}
static inline int print_dnarc(uint64_t dna, uint64_t rc, bool dbg = true)
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
static inline uint64_t revcmp(uint64_t dna)
{
    uint64_t m = 0x3333333333333333UL;
    dna = ((dna & m) << 2) | ((dna >> 2) & m);
    m = 0x0f0f0f0f0f0f0f0fUL;
    dna = ((dna & m) << 4) | ((dna >> 4) & m);
    asm ("bswap %0" : "=r" (dna) : "0" (dna));
    dna ^= 0xaaaaaaaaaaaaaaaaUL;
    return dna >> (64 - (KEY_WIDTH << 1));
}

static void show_list(kct_t* kc, std::list<uint32_t> &bnd)
{
    std::list<uint32_t>::iterator it = bnd.begin();
    unsigned i = 0;
    for (it = bnd.begin(); it != bnd.end(); ++it) {
        Bnd* bd = &kc->bd[*it];
        EPR("[%u]:\t%u\t%u", i++, bd->s, bd->l);
    }
}

static inline int print_ndx(uint64_t dna, bool dbg = true)
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
enum ensembl_parts {ID,  SEQTYPE, IDTYPE,
        IDTYPE2, BUILD, ID2, START, END, NR, META};

enum bnd_type {NEW_REF_ID, START_SEQ, N_STRETCH, UQ_REGION, END_REF, DBSNP, KNOWN_SITE};

//ensembl format: >ID SEQTYPE:IDTYPE LOCATION [META]
// fai does not handle chromosomes with offset.
Hdr* new_header(kct_t* kc, void* g, int (*gc) (void*))
{
    int c;
    Hdr* h = new Hdr, dbg;
    kc->h.push_back(h);
    h->hdr_type = ENSEMBL_HDR;

    uint32_t p = ID;
    h->part[p] = kc->id_l;
    while ((c = gc(g)) != -1) {
        //fputc(c, stderr);
        switch (c) {
            case '\n': _buf_grow_add_err(kc->id, 1ul, '\0', return NULL); break;
            case ' ':
                _buf_grow_add_err(kc->id, 1ul, '\0', return NULL);
                if (p == IDTYPE) {
                    char* q = &kc->id[h->part[p]];
                    if (strncmp(q, "chromosome", strlen(q)) != 0 &&
                            strncmp(q, "supercontig", strlen(q)) != 0 &&
                            strncmp(q, "nonchromosomal", strlen(q)) != 0)
                        h->hdr_type = UNKNOWN_HDR;
                } else if (p != ID) {
                    h->hdr_type = UNKNOWN_HDR;
                }
                if (p < ARRAY_SIZE(h->part))
                    h->part[++p] = kc->id_l;
                continue;
            case ':': 
                _buf_grow_add_err(kc->id, 1ul, '\0', return NULL);
                if (p == ID || p == IDTYPE)
                    h->hdr_type = UNKNOWN_HDR;
                else if (p == SEQTYPE) {
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
                        } else { // could still be ok.
                            h->hdr_type = UNKNOWN_HDR;
                        }
                    }
                }
                if (p < ARRAY_SIZE(h->part))
                    h->part[++p] = kc->id_l;
                continue;
            default:
                if (isspace(c)) {
                    c = '\0';
                    h->hdr_type = UNKNOWN_HDR;
                }
                _buf_grow_add_err(kc->id, 1ul, c, return NULL);
                continue;
        }
        break;
    }
    if (p < NR) 
        h->hdr_type = UNKNOWN_HDR;

    kc->hdr.insert(std::pair<char*, Hdr*>(kc->id + h->part[ID], h));
    _buf_grow_err(kc->bd, 1ul, return NULL);
    ASSIGN_BD(kc->bd[kc->bd_l], NEW_REF_ID, ~0u, 0, ~0u, 0, 0);

    if (h->hdr_type == ENSEMBL_HDR) {
        h->p_l = p;
        p = atoi(kc->id + h->part[START]);
        if (p != 1)
            kc->bd[kc->bd_l].s = p;
        kc->bd[kc->bd_l].l = atoi(kc->id + h->part[END]);
    } else {
        h->p_l = 1;
        EPR("\nWARNING: non-ensembl reference.");
    }
    return h;
}

void free_kc(kct_t* kc)
{
    for (Kct* k = kc->kct; k != &kc->kct[kc->kct_l]; ++k)
        if (k->seq.m > 3)
            free(k->p.b2);

    std::map<char*, Hdr*, Tid>::iterator it;
    for (it = kc->hdr.begin(); it != kc->hdr.end(); ++it)
        delete it->second;

    _buf_free(kc->id);
    _buf_free(kc->bd);
    _buf_free(kc->kct);
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
    uint32_t b2pos = 0u, addpos = 0u, endpos = 0u;
    int c, dbg = 0;
    uint32_t t = KEY_WIDTH;
    Kct *ct = NULL;
    uint8_t *q = NULL;
    Hdr* h = NULL;
    const char* hdr = 0;

    while ((c = gc(g)) != '>' && c != '@' && c != -1) {} // skip to first ref ID

    // TODO: realpos based dbsnp or known site boundary insertion.
    while (c >= 0) {
        if (c != '>') { // N-stretch
            t = 0; // use  to count stretch (key is rebuilt later)
            do {
                if (c != '\n') {
                    if (c != 'N') {
                        if (B6(b, c) || c == '>') break;
                        EPR("Ignoring strange nucleotide:%c\n", c);
                    }
                    ++t;
                }
                c = gc(g);
            } while(c != -1);
            if (c == -1) break;
            ungc(c, g);
            _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
            ASSIGN_BD(kc->bd[kc->bd_l], N_STRETCH, ~0u, b2pos + addpos, t + KEY_WIDTH, 0, dna);
            addpos += t;
        } else { // header
            if (h != NULL) {
                // add enddna as length to first element, the chromome boundary.
                DEBUG_ASSIGN_ENDDNA(kc->bd[*(h->bnd.begin())].end_dna, dna);
                _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
                ASSIGN_BD(kc->bd[kc->bd_l], END_REF, ~0u, b2pos + addpos, 0, 0, dna);
                h->bnd.push_back(kc->bd_l++);
            }
            b2pos = 0u;
            h = new_header(kc, g, gc);
            hdr = kc->id + h->part[0];
            if (h == NULL) return -EFAULT;
            addpos = kc->bd[kc->bd_l].s;
            endpos = kc->bd[kc->bd_l].l;
        }
        for (t = b2pos + KEY_WIDTH; b2pos != t; ++b2pos) {
            b ^= b;
            switch(c = gc(g)) {
                case 'C': case 'c': b = 2;
                case 'G': case 'g': b ^= 1;
                case 'U': case 'u': case 'T': case 't': b ^= 2;
                case 'A': case 'a': dna = _seq_next(b, dna, rc);
                case '\n': continue;
            }
            break;
        }
        if_ever (b2pos != t) {// c == -1, 'N' or '>'
            EPR("key incomplete before next boundary - "
                    "may cause enddna assertion failure");
            continue;
        }
        // ndx (key) is 2nd cNt bit excised, rc or dna
        _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
        kc->bd[kc->bd_l].dna = dna;
        //EPR("dna:0x%lx\trc:0x%lx", dna, rc);
        //print_dna(dna);
        h->bnd.push_back(kc->bd_l);
        do {
            ndx = _get_ndx(ndx, dna, rc);
            ASSERT(ndx < (1ul << KEYNT_BUFSZ_SHFT), return -EINVAL)
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
            } else { //EPR("new key 0x%lx at %u", ndx, b2pos);
                _buf_grow(kc->kct, 1ul);
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
                    //EPQ(ndx == 0x47ef9, "%c added to (%x) <===== ndx:0x%lx dna:0x%lx, %u", b6(b << 1), *q, ndx, dna, b2pos);
                    dna = _seq_next(b, dna, rc);
                    /*if ((strncmp("GL000247.1", hdr, strlen(hdr)) == 0 && b2pos > 5200 && b2pos < 5300) ||
                        (strncmp("GL000244.1", hdr, strlen(hdr)) == 0 && b2pos > 38600 && b2pos < 38700) ||
                        (strncmp("GL000243.1", hdr, strlen(hdr)) == 0 && b2pos > 12000 && b2pos < 12100) ||
                        (strncmp("GL000194.1", hdr, strlen(hdr)) == 0 && b2pos > 1750 && b2pos < 1850)) {
                            EPR0("%s\t%u\t", hdr, b2pos);
                            print_dna(dna);
                    }*/
                    ++b2pos;
                    continue;
            }
            if (ct->seq.m == 3) --ct->seq.l;
            else --ct->p.l;
            break;
        } while (c >= 0);

        if (c == '>' || c == -1) {
            EPR("processed %u Nts for %s", b2pos, kc->id + h->part[0]);
            if (b2pos + addpos != endpos) {
                t = *(h->bnd.begin());
                EPR("b2pos + addpos != endpos: %u + %u == %u (%u)", b2pos, addpos, endpos, kc->bd[t].l);
                print_dna(dna);
                // correct it.
                //endpos = b2pos + addpos;
                kc->bd[t].l = b2pos + addpos;
            }
        } else {
            EPR("=>\tN-stretch at Nt %u", b2pos);
        }
        DEBUG_ASSIGN_ENDDNA(kc->bd[kc->bd_l].end_dna, dna);
        ++kc->bd_l;
    }
    if (h == NULL)
        return -EFAULT;
    // add enddna as length to first element, the chromome boundary.
    DEBUG_ASSIGN_ENDDNA(kc->bd[*(h->bnd.begin())].end_dna, dna);
    _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
    ASSIGN_BD(kc->bd[kc->bd_l], END_REF, ~0u, b2pos + addpos, 0, 0, dna);
    h->bnd.push_back(kc->bd_l++);
    //EPR("%s, %u stored", str, bd.len & NR_MASK);
    return 0;
}

int decr_excise(kct_t* kc, uint32_t* wbuf, Walker* wlkr, unsigned i, unsigned left)
{
    while (--i != left) {
        //EPR0("%u/%u ", i, left);
        ASSERT(wbuf[i] != ~0u, return -EFAULT, "%u/%u", i, left);
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
            ASSERT(((x->seq.l-1ul) >> 2) < ARRAY_SIZE(x->seq.b2),
                  return -EFAULT);
            qe = &x->seq.b2[(x->seq.l-1) >> 2];
            ASSERT(x->seq.l != 0, return -EINVAL);
            t = --x->seq.l;
        } else {
            q = &x->p.b2[t >> 2];
            qe = &x->p.b2[(x->p.l-1) >> 2];
            // could convert to Kct.seq, but why bother?
            ASSERT(x->p.l != 0, return -EINVAL);
            t = --x->p.l;
        }
        if (t == 0) {
            EPR0("two successive decrements on same key:\t");
            print_ndx(wbuf[i]);
        }
        wbuf[i] = ~0u;
        //EPQ0(dbg, "%c", b6(((*q >> ((t & 3) << 1)) << 1) & 6));
        //EPR(dbg, "decreasing:%x at %lu:%c", (unsigned)*q, t,
        //        b6(((*q >> ((t & 3) << 1)) << 1) & 6));
        t = (1 << (((w->count & 3) + 1) << 1)) - 1;
        *q = ((*q & (t ^ 0xff)) >> 2) | (*q & (t >> 2));
        //EPR("became:%x", (unsigned)*q);
        while (q != qe) {
            //EPR("loop start:%x", (unsigned)*q);
            *q |= (q[1] & 3) << 6;
            //EPR("loop end:%x", (unsigned)*q);
            *++q >>= 2;
        }
        //EPR("and last:%x", (unsigned)*q);
    }
#ifdef DEBUG
    wbuf[i] = ~0u;
    while (i != 0) {
        if (wbuf[i] != ~0u) {
            Walker* w = &get_w(wlkr, kc, wbuf[i]);
            ASSERT(wbuf[i] == ~0u, return print_ndx(wbuf[i]), "index %u/%u:%u", i, left, w->tmp_count);
        }
        --i;
    }
#endif
    return left;
}

// XXX: inferiority ok? XXX: store deviant bit when unique,
// reverse seq could be required for mapping - maybe 
int extend_uniq(kct_t* kc, const int ext)
{
    uint64_t t = kc->kct_l * sizeof(Walker), dna = 0ul;
    Walker* wlkr = (Walker*)malloc(t);
    memset(wlkr, 0ul, t);
    uint32_t uqct, iter = 0;
    uint32_t* wbuf = (uint32_t*)malloc(ext * sizeof(uint32_t));
    for (int i = 0; i != ext; ++i) wbuf[i] = ~0u;
    bool dbg = false;
    const char* dbgtid = "GL000246.1 ";
    // make sure we have always one space available;
    _buf_grow_err(kc->bd, 1ul, return -ENOMEM);

    do { // until no no more new uniques
        uqct = 0;
        std::list<Hdr*>::iterator h;
        for (h = kc->h.begin(); h != kc->h.end(); ++h) { // over ref headers

            std::list<uint32_t>::iterator bdit = (*h)->bnd.begin();
            ASSERT(bdit != (*h)->bnd.end(), return -EFAULT);

            // first is special: chromo, which marks start and end of chromosome
            Bnd* bd = &kc->bd[*bdit];
            uint32_t pos = bd->s, endpos = bd->l;
            const char* hdr = kc->id + (*h)->part[0];
            EPR("----[\t%s:%u-%u\t]----", hdr, pos, endpos);
            pos += KEY_WIDTH;
            while (++bdit != (*h)->bnd.end()) {
                if (strncmp(hdr, dbgtid, strlen(dbgtid)) == 0) {
                    EPR("... turning debugging on ...");
                    dbg = true;
                }// else { dbg = false; }
                // loop over events (stretches and maybe later, SNVs, splice sites)

                dna = bd->dna;
                uint64_t rc = revcmp(dna);
                uint64_t ndx = _getxtdndx(kc, ndx, dna, rc);
                uint32_t bd_i = kc->bd_l;
                uint32_t infior = 1;
                EPR("dna:0x%lx\trc:0x%lx,%u/%u", dna, rc, *bdit, bd_i);
                print_dna(dna);
                kc->bd[bd_i].s = pos;
                kc->bd[bd_i].t = START_SEQ;
                kc->bd[bd_i].end_dna = dna;
                bd = &kc->bd[*bdit];

                //show_list(kc, (*h)->bnd);

                ASSERT(pos - KEY_WIDTH < endpos, return -EINVAL, "%u < %u", pos, endpos);

                int left = ext;
                EPR("Next event at %u, now %u", bd->s, pos);
                ASSERT(pos <= bd->s, return -EFAULT);

                for (;pos != bd->s;++pos) { // until next event
                    ASSERT(_getxtdndx0(kc, ndx) != UNINITIALIZED, return -print_dna(dna), "at %u, 0x%lx", pos, ndx);

                    Kct *y = &_get_kct(kc, ndx);
                    Walker* w = &get_w(wlkr, kc, ndx);
                    if (left == ext) {
                        //EPR("position after unique, unless extended.");
                        // what we'll jump to, in subsequent iterations.
                        // This is only inserted when bd_i increments -
                        // otherwise it's overwritten successively.
                        kc->bd[bd_i].dna = dna;
                        kc->bd[bd_i].l = pos - kc->bd[bd_i].s;
                        wbuf[--left] = ndx;
                        if (infior > w->infior)
                            w->infior = infior;

                    } else if (left) { // keep filling the buffer
                        if (--left) {
                            wbuf[left] = ndx;
                            if (infior > w->infior)
                                w->infior = infior;
                        } else { // actually buf just became full, without extension or already handled..
                            infior = 0;

                            if (bd_i != kc->bd_l) {
                                ASSERT(bd_i + 1 == kc->bd_l, return -1);
                                // we did insert a range
                                EPQ(dbg, "range ended at\t%s\t%u", kc->id + (*h)->part[0],
                                        kc->bd[bd_i].s + kc->bd[bd_i].l);
                                print_dna(dna, dbg);
                                //EPR("before %u, bd_I:%u", *bdit, bd_i);
                                //show_list(kc, (*h)->bnd);
                                ASSERT(bdit != (*h)->bnd.end(), return -EFAULT);
                                (*h)->bnd.insert(bdit, bd_i++);
                                //EPR("after %u", *bdit);
                                //show_list(kc, (*h)->bnd);
                                _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
                                bd = &kc->bd[*bdit];
                            }
                            // add tmp_count for each - we cannot extend it in this iteration
                            for (left = ext; --left;) {
                                get_w(wlkr, kc, wbuf[left]).count++;
                                ASSERT(get_w(wlkr, kc, wbuf[left]).tmp_count > 0, return -EFAULT);
                                --get_w(wlkr, kc, wbuf[left]).tmp_count;
                                wbuf[left] = ~0u;
                            }
                            ASSERT(w->tmp_count == 0, return -EFAULT, "(%u)", w->tmp_count);

                        }
                    }
                    t = w->count + w->tmp_count;
                    if (y->seq.m == 3) { //EPR("Kct in seq format");
                        ASSERT(t < y->seq.l, return -EINVAL, "%s\t%u", kc->id + (*h)->part[0], pos);
                        t = (y->seq.b2[t >> 2] >> ((t & 3) << 1)) & 3;
                        dna = _seq_next(t, dna, rc);
                        EPQ0(dbg, "ndx:0x%lx at %u. dna|rc:\t", ndx, pos);
                        print_dnarc(dna, rc, dbg);
                        t = y->seq.l;
                    } else {
                        //EPQ0(dbg, "Kct in p format. m:%lu, l:%lu, t:%lu, ndx:0x%lx\t0x%lx\t", y->p.m, y->p.l, t, ndx, (size_t)&y->p.b2); print_dna(dna, dbg);
                        ASSERT(t < y->p.l, return -EINVAL);

                        t = (y->p.b2[t >> 2] >> ((t & 3) << 1)) & 3; // get next 2bit

                        dna = _seq_next(t, dna, rc);
                        EPQ0(dbg, "dna:\t"); print_dna(dna, dbg);

                        //fprintf(stderr, "%c\n", b6(t<<1));
                        t = y->p.l; // can also be 1 (unique), unless we decide to convert back upon decrement
                        //EPR("end Kct in p format. m:%lu, l:%lu, t:%lu, ndx:0x%lx\t0x%lx", y->p.m, y->p.l, t, ndx, (size_t)&y->p.b2);
                    }
                    ndx = _getxtdndx(kc, ndx, dna, rc);
                    if (t > 1ul) {
                        if (left == 0) ++w->count;
                        else ++w->tmp_count;
                        continue;
                    }
                    ASSERT(t == 1, return -EINVAL);
                    // found unique
                    if (left == 0) { // this is a first unique.
                        EPQ(dbg, "unique at\t%s\t%u", kc->id + (*h)->part[0], pos);
                        print_dna(dna, dbg);
                        ASSIGN_BD(kc->bd[kc->bd_l], UQ_REGION, ~0u, pos, ~0u, infior, dna);
                        left  = ext;
                        infior = w->infior + 1;
                        continue;
                    }
                    if (bd_i == kc->bd_l) { // a second unique just now.
                        // insert a new range.
                        EPQ(dbg, "new range till\t%s\t%u", kc->id + (*h)->part[0], pos);
                        ++uqct;
                        ++kc->bd_l; // incr. to invoke insert when range completes.
                        left = decr_excise(kc, wbuf, wlkr, ext, left);
                        if (left < 0)
                            return left;
                        left = ext;
                        continue;
                    }
                    // else another unique within range.
                    EPQ(dbg, "range extended to\t%s\t%u", kc->id + (*h)->part[0], pos);

                    // XXX: infior can be more conservatively increased once an
                    // original range boundary is reached. Extended with another
                    // range.
                    ++infior;

                    //
                    ASSERT(left != ext, return -EFAULT);
                    left = decr_excise(kc, wbuf, wlkr, ext, left);
                    if (left < 0)
                        return left;
                    //w = &get_w(wlkr, kc, wbuf[left]);
                    //ASSERT(w->tmp_count == 0, return -EFAULT);
                    EPQ0(dbg, "\n");
                    print_dna(dna, dbg);
                    left = ext;

                    // XXX: should inferiority increase?
                    // found second unique
                }
                if (left && left != ext) {
                    if (bd_i != kc->bd_l) {
                        EPR("bd_i(%u) != kc->bd_l(%u)", bd_i, kc->bd_l);
                        // we did insert a range
                        //(*h)->bnd.insert(bdit, bd_i++);
                        //show_list(kc, (*h)->bnd);
                        //_buf_grow_err(kc->bd, 1ul, return -ENOMEM);
                        //bd = &kc->bd[*bdit];
                    }
                    // add tmp_count for each - we cannot extend it in this iteration
                    /*for (left = ext; --left;) {
                        get_w(wlkr, kc, wbuf[left]).count++;
                        ASSERT(get_w(wlkr, kc, wbuf[left]).tmp_count > 0, return -EFAULT);
                        --get_w(wlkr, kc, wbuf[left]).tmp_count;
                        wbuf[left] = ~0u;
                    }*/
                    EPR("range extended until boundary");
                    left = decr_excise(kc, wbuf, wlkr, ext, --left);
                    if (left < 0)
                        return left;
                }
#ifdef DEBUG
                if (kc->bd[*bdit].t != N_STRETCH) {
                    ASSERT(kc->bd[*bdit].end_dna == dna,
                            return -print_dnarc(dna, rc, 1) && print_dna(kc->bd[*bdit].end_dna));
                }
#endif
                EPR("%u => %u", pos, pos + kc->bd[*bdit].l);
                pos += kc->bd[*bdit].l; // add skipped Nts
                //if (kc->bd[*bdit].t == N_STRETCH) {
                //    ++pos;
                //}
                EPR0("end of range: ndx:0x%lx\tdna|rc:(%u)\t", ndx, pos);
                print_dnarc(dna, rc, 1);
            }
        }
        EPR("extended %u unique ranges in iteration %u", uqct, ++iter);
        for (t = 0; t != kc->kct_l; ++t) {
            Walker* w = wlkr + t;
            ASSERT(w->tmp_count == 0u, return -EFAULT, "%lu", t);
            w->count = 0u;
        }
        dbg = true;
    } while (uqct != 0);
    free(wlkr);
    free(wbuf);
    return 0;
}

int
fa_print(seqb2_t *fa)
{
    return 1;
}

int
fn_convert(struct gzfh_t* fhout, const char* search, const char* replace)
{
    char* f = strstr(fhout->name, search);
    if (f == NULL) return 0;
    strncpy(f, replace, strlen(replace) + 1);
    return 1;
}

int fa_index(struct seqb2_t* seq)
{
    struct gzfh_t* fhout = seq->fh + ARRAY_SIZE(seq->fh) - 1;
    struct gzfh_t* fhin = seq->fh + 2; // init with ref
    uint64_t blocksize = (uint64_t)seq->blocksize << 20;
    int res, ret = -1;
    void* g;
    int (*gc) (void*);
    int (*ungc) (int, void*);
    int is_gzfile = fhin->io != NULL;
    char file[256];

    const char* ndxact[4] = {".fa", "_kcpos.wig.gz"};

    if (fhout->name != NULL && ((res = strlen(fhout->name)) > 255)) {
        strncpy(file, fhout->name, res);
    } else {
        res = strlen(fhin->name);
        ASSERT(strncmp(fhin->name, "stdin", res) != 0, return -1);
        strcpy(file, fhin->name);
    }
    fhout->name = &file[0];

    kct_t kc;
    kc.kcsndx = _buf_init_arr(kc.kcsndx, KEYNT_BUFSZ_SHFT);
    kc.kct = _buf_init(kc.kct, 16);
    kc.bd = _buf_init(kc.bd, 0);

    // FIXME: assigning 1 Mb for reference id, to keep realloc from happening.
    // (realloc would invalidate pointers inserted in kc->hdr)
    kc.id = _buf_init(kc.id, 20);

    // the first entry (0) is for a reference ID always.
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
    res = extend_uniq(&kc, seq->readlength - KEY_WIDTH);
    ASSERT(res >= 0, goto out);

    ret = res != 0;

out:
    free_kc(&kc);
    return ret;

}


