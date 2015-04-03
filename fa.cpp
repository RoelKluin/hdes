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
        EPR0("[%u]:\t%u\t%u\t", i++, bd->s, bd->l);
        print_dnarc(bd->at_dna, bd->dna);
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

enum bnd_type {NEW_REF_ID, N_STRETCH, UQ_REGION, END_REF, DBSNP, KNOWN_SITE};

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
    // NB: bij NEW_REF_ID zijn ->s en ->l andersom.
    ASSIGN_BD(kc->bd[kc->bd_l], NEW_REF_ID, UNINITIALIZED, UNINITIALIZED, 0, 0, ~0u);

    if (h->hdr_type == ENSEMBL_HDR) {
        h->p_l = p;
        // NB: bij NEW_REF_ID zijn ->s en ->l andersom.
        kc->bd[kc->bd_l].l = atoi(kc->id + h->part[START]) - 1;
        kc->bd[kc->bd_l].s = atoi(kc->id + h->part[END]);
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
            ASSIGN_BD(kc->bd[kc->bd_l], N_STRETCH, UNINITIALIZED, b2pos + addpos, t, 0, dna);
            addpos += t;
        } else { // header
            if (h != NULL) {
                // add enddna as length to first element, the chromome boundary.
                //DEBUG_ASSIGN_ENDDNA(kc->bd[*(h->bnd.begin())].at_dna, dna);
                _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
                ASSIGN_BD(kc->bd[kc->bd_l], END_REF, UNINITIALIZED, b2pos + addpos, 0, 0, dna);
                h->bnd.push_back(kc->bd_l++);
            }
            b2pos = 0u;
            h = new_header(kc, g, gc);
            hdr = kc->id + h->part[0];
            if (h == NULL) return -EFAULT;
            // NB: bij NEW_REF_ID zijn ->s en ->l andersom.
            addpos = kc->bd[kc->bd_l].l;
            endpos = kc->bd[kc->bd_l].s;
        }
        for (t = b2pos + KEY_WIDTH; b2pos != t; ++b2pos) {
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
        // ndx (key) is 2nd cNt bit excised, rc or dna
        _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
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
                    dna = _seq_next(b, dna, rc);
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
                // correct it.
                // NB: bij NEW_REF_ID zijn ->s en ->l andersom.
                kc->bd[t].s = b2pos + addpos;
            }
        } else {
            EPR("=>\tN-stretch at Nt %u", b2pos);
        }
        //DEBUG_ASSIGN_ENDDNA(kc->bd[kc->bd_l].at_dna, dna);
        ++kc->bd_l;
    }
    if (h == NULL)
        return -EFAULT;
    // add enddna as length to first element, the chromome boundary.
    //DEBUG_ASSIGN_ENDDNA(kc->bd[*(h->bnd.begin())].at_dna, dna);
    _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
    ASSIGN_BD(kc->bd[kc->bd_l], END_REF, UNINITIALIZED, b2pos + addpos, 0, 0, dna);
    h->bnd.push_back(kc->bd_l++);
    return 0;
}

int decr_excise(kct_t* kc, uint32_t* wbuf, Walker* wlkr, unsigned i, unsigned left)
{
    while (i-- != left) {
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
            t = --x->seq.l;
        } else {
            q = &x->p.b2[t >> 2];
            qe = &x->p.b2[(x->p.l-1) >> 2];
            // could convert to Kct.seq, but why bother?
            t = --x->p.l;
        }
        t = (1 << (((w->count & 3) + 1) << 1)) - 1;
        *q = ((*q & (t ^ 0xff)) >> 2) | (*q & (t >> 2));
        wbuf[i] = UNINITIALIZED;
        while (q != qe) {
            *q |= (q[1] & 3) << 6;
            *++q >>= 2;
        }
    }
    wbuf[i] = UNINITIALIZED;
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
    for (int i = 0; i != ext; ++i) wbuf[i] = UNINITIALIZED;
    int dbg = 0;
    const char* dbgtid = "MT ";//GL000239.1";//GL000210.1"; //GL000239.1";//MT ";
    // make sure we have always one space available;
    _buf_grow_err(kc->bd, 1ul, return -ENOMEM);

    do { // until no no more new uniques
        uqct = 0;
        std::list<Hdr*>::iterator h;
        for (h = kc->h.begin(); h != kc->h.end(); ++h) { // over ref headers

            std::list<uint32_t>::iterator bdit = (*h)->bnd.begin();

            // first is special: chromo, which contains start and end of chromosome
            // l en s zijn andersom opdat bij extensie de juiste, l, groeit
            Bnd *last, *next = &kc->bd[*bdit];
            uint32_t pos = 0, add = 0;//, endpos = next->s;
            const char* hdr = kc->id + (*h)->part[0];

            while (++bdit != (*h)->bnd.end()) {
                // N-stretches, later also skips, future: SNVs, splice sites?
                last = next;
                next = &kc->bd[*bdit];
                dna = last->dna; // first seq after skip
                pos += last->l; // gaat ook goed voor eerste 'last', van ref.
                dbg = strncmp(hdr, dbgtid, strlen(dbgtid)) == 0;

                //top = last; // XXX


                uint64_t rc = revcmp(dna);
                uint64_t ndx = _getxtdndx(kc, ndx, dna, rc);
                if (last->t == N_STRETCH) {
                    pos += KEY_WIDTH;
                    //add += last->l;
                } else if (last->t == NEW_REF_ID) {
                    pos += KEY_WIDTH;
                } else {
                    EPR("type: %u", last->t);
                }

                Bnd *top = &kc->bd[kc->bd_l];
                top->l = UNINITIALIZED;
                top->s = pos;
                top->at_dna = ~0u;

                uint32_t infior = 1;
                // boundary is considered as a first unique,
                // if we find a 2nd, we really want to extend this region
                int left = ext;
                EPR0("----[\t%s%s:%u-%u\t\t]----\tndx:", strlen(hdr) < 8 ? "\t":"",hdr, pos, next->s); print_ndx(ndx);

                for (;pos < next->s + add;++pos) { // until next event
                    ASSERT(_getxtdndx0(kc, ndx) != UNINITIALIZED, return -print_dna(dna), "at %u, 0x%lx", pos, ndx);
                    print_dna(dna, dbg);

                    Kct *y = &_get_kct(kc, ndx);
                    Walker* w = &get_w(wlkr, kc, ndx);
                    if (left == ext) {
                        //EPR("position after 2nd unique, or each extended.");
                        // what we'll jump to, in subsequent iterations.
                        // This is only inserted when bd_i increments -
                        // otherwise it's overwritten successively.
                        top->dna = dna;
                        top->l = pos - top->s;
                        wbuf[--left] = ndx;
                        if (infior > w->infior)
                            w->infior = infior;

                    } else if (left) { // keep filling the buffer
                        if (--left) {
                            wbuf[left] = ndx;
                            if (infior > w->infior)
                                w->infior = infior;
                        } else { // actually buf just became full, without extension(s not already handled)..
                            EPQ(top->s > 0xfffffff, "%u", -top->s);
                            infior = 0;

                            if (top->l != UNINITIALIZED) {
                                // we did insert a range
                                print_dnarc(top->at_dna, top->dna, dbg);
                                //show_list(kc, (*h)->bnd);
                                (*h)->bnd.insert(bdit, kc->bd_l++); // evan after increment one should be available for top.
                                _buf_grow_err(kc->bd, 1ul, return -ENOMEM);
                                last = NULL;
                                top = &kc->bd[kc->bd_l];
                                top->l = UNINITIALIZED;
                                next = &kc->bd[*bdit]; // unchanged.
                            }
                            // add tmp_count for each - we cannot extend it in this iteration
                            for (left = ext; --left;) {
                                get_w(wlkr, kc, wbuf[left]).count++;
                                --get_w(wlkr, kc, wbuf[left]).tmp_count;
                                wbuf[left] = UNINITIALIZED;
                            }
                        }
                    }
                    t = w->count + w->tmp_count;
                    if (y->seq.m == 3) { //EPR("Kct in seq format");
                        ndx = (y->seq.b2[t >> 2] >> ((t & 3) << 1)) & 3;
                        t = y->seq.l;
                    } else {
                        ndx = (y->p.b2[t >> 2] >> ((t & 3) << 1)) & 3; // get next 2bit
                        t = y->p.l; // can also be 1 (unique), unless we decide to convert back upon decrement
                    }
                    dna = _seq_next(ndx, dna, rc);
                    ndx = _getxtdndx(kc, ndx, dna, rc);
                    if (t > 1ul) {
                        if (left == 0) ++w->count;
                        else ++w->tmp_count;
                        continue;
                    }
                    ++w->tmp_count;
                    // found unique
                    ++infior;
                    if (left == 0) { // this is a first unique.
                        top = &kc->bd[kc->bd_l];
                        top->t = UQ_REGION;
                        top->s = pos;
                        top->i = max(infior, w->infior + 1u);
                        top->at_dna = dna;
                        left = ext;
                        infior = top->i + 1;
                        continue;
                    }
                    // insert a new range.
                    ++uqct;
                    left = decr_excise(kc, wbuf, wlkr, ext, left);
                    if (left < 0) return left;
                    left = ext;
                }
                EPQ(dbg, "---");
                ASSERT(dna == next->at_dna, return -print_dnarc(dna, next->at_dna), "=> enddna mismatch");
                if (last && left && left != ext) {
                    EPR("joining region");

                    last->l += last->s - next->s;
                    bdit = (*h)->bnd.erase(bdit);
                    --kc->bd_l;
                    top = &kc->bd[kc->bd_l];
                    top->s = -3u;
                    // add tmp_count for each - we cannot extend it in this iteration
                    left = decr_excise(kc, wbuf, wlkr, ext, left);
                    if (left < 0)
                        return left;
                }
//                pos += kc->bd[*bdit].l; // add skipped Nts
                EPQ0(dbg, "%u => %u ]\t\t\t\tndx:\t", kc->bd[*bdit].s, kc->bd[*bdit].s + kc->bd[*bdit].l);print_ndx(ndx,dbg);
            }
            show_list(kc, (*h)->bnd);
        }
        //dbg = 1;
        EPR("extended %u unique ranges in iteration %u", uqct, ++iter);
        for (t = 0; t != kc->kct_l; ++t) {
            Walker* w = wlkr + t;
            w->count = 0u;
        }
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


