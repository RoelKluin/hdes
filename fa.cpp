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
#include <limits.h> //INT_MAX
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
//#include <glib.h>
//#include <pcre.h>
#include "fa.h"

static int
start_dbg(kct_t *C kc, Hdr *C h, uint64_t C dna)
{
    char C* hdr = kc->id + h->part[0];

    // TODO: if not start of contig or just after a N-stretch, the dna key should be
    // just after an unique key. we could check this.

    dbg = strncmp(hdr, dbgchr, strlen(hdr)) ? 3 : 5;

    if (dbg > 3) {
        EPR0("----[\t%s%s:(%u+)%u..%u\t]----\tdna:", strlen(hdr) < 8 ? "\t": "",
                hdr, (*kc->bdit).corr, (*kc->bdit).s, (*kc->bdit).e);
        print_dna(dna);
        show_mantras(kc, h);
    }
    return 0;
}

void
free_kc(kct_t* kc)
{
    std::list<Hdr*>::iterator it;
    for (it = kc->h.begin(); it != kc->h.end(); ++it) {
        Hdr* h = *it;
        free(h->part);
        delete h;
    }
    _buf_free(kc->id);
    _buf_free(kc->s);
    _buf_free(kc->kct);
    _buf_free(kc->ndxkct);
}

static inline void
update_max_infior(uint64_t *C k, uint64_t infior)
{
    expect(infior < MAX_INFIOR) {
        infior += INFIOR;
    } else {
        WARN("MAX_INFIOR reached");
    }
    if (infior > *k)
        *k ^= (*k ^ infior) & INFIOR_MASK;
}

static inline C int
get_nextnt(kct_t C*C kc, uint64_t p)
{
    ASSERT(p < (kc->s_l << 2), return -EFAULT, "%lu/%lu", p, kc->s_l);
    return (kc->s[p>>2] >> ((p&3) << 1)) & 3;
}

static inline int
mark_leaving(kct_t C*C kc, uint64_t *C k, uint64_t p)
{
    int res = -EFAULT;
    uint64_t nt;
    if (IS_DBG_K(kc, k))
        DESCRIBE_KEY(kc, k, '<');
    if (p == ~0ul) { //end of contig or before N-stretch or boundary: no nextNt.
        if (IS_DISTINCT(k)) {
            EPQ(IS_DBG_K(kc, k), " ==> distinct cleared");
            *k &= ~(DISTINCT|NT_MASK); // clear from last iteration
        }
        EPQ(IS_DBG_K(kc, k), "key at end of contig, before uq-bnd or N-stretch always marked as distinct");
        *k |= DISTINCT;
    } else if (ALL_SAME_NTS(k) == false) { // formerly PENDING(k), this should be ok?

        nt = (kc->s[p>>2] >> ((p&3) << 1)) & 3;

        EPQ(IS_DBG_K(kc, k), "\t%c", b6(nt<<1));

        if (IS_FIRST(k)) {
            if (IS_DISTINCT(k)) {
                EPQ(IS_DBG_K(kc, k), " ==> distinct cleared");
                *k &= ~(DISTINCT|NT_MASK); // clear from last iteration
            } else if (kc->iter) {
                //if ok, then we can remove the kct iteration post loop?
//                ASSERT(*k & MARKED, goto err, "MARKED not obsolete"); // TODO: check after setting marked
            }
            *k |= nt << NT_SHFT;
        } else {
            uint64_t first_nt = FIRST_NT(k);
            EPQ(IS_DBG_K(kc, k), "%lu:%c <=> %c", K_OFFS(kc, k), b6(first_nt<<1), b6(nt<<1)); 
            if (IS_DISTINCT(k) == false && nt != first_nt) {
                EPQ(IS_DBG_K(kc, k), " ==> dbgtsoffs got distinct");
                *k |= DISTINCT;
            }
        }
    }

    if (IS_LAST(k)) {
	*k &= ~AT_MASK; /* reset position for next iterationtracking */
        *k ^= DISTINCT;
    }
    res = 0;
err:
    if (res || IS_DBG_K(kc, k))
	DESCRIBE_KEY(kc, k, '?');

    return res;
}

static int
_mark_all(kct_t C*C kc, char C*C f, unsigned C l, unsigned C len, unsigned C end)
{
    unsigned p = end - len, i = 0;
    int res = 0;
    while (p < end && res == 0) {
        uint64_t* k = kc->kct_scope[i++];
        EPQ(IS_DBG_K(kc, k), "%s +%u mark_leaving()", f, l);
        res = mark_leaving(kc, k, p++);
    }
    if (res < 0)
        EPQ(res, "error during mark_all():%s +%u", f, l);
    return res;
}
#define mark_all(kc, i, e) _mark_all(kc, __FILE__, __LINE__, i, e)

static inline void
decr_excise(kct_t *C kc, uint64_t C*C kct, uint64_t C infior, unsigned C i)
{
    uint64_t *C k = kc->kct_scope[i];
    if (k != kct) // .. don't elevate key if it is the one that became uniq
        update_max_infior(k, infior);
    EPQ(IS_DBG_K(kc, k), "exiced: [%u]", K_OFFS(kc, k));
    ++kc->uqct;
    if ((--*k & UQ_MASK) == 1ul) {

        *k ^= UQ_BIT | 1;
 
    } else if (IS_LAST(k)) {
        *k &= ~AT_MASK;
    }
}

static void
decr_excise_all(kct_t *C kc, uint64_t * k, uint64_t C*C kct, unsigned const end)
{
    if (k == NULL) return;
    uint64_t infior = (kct ? std::max(*k, *kct) : *k) & MAX_INFIOR;
    expect(infior < MAX_INFIOR) {
        infior += INFIOR;
    } else {
        WARN("MAX_INFIOR reached");
    }
    unsigned i = IS_UQ(k) ? 1 : 0;
    while (i < end)
        decr_excise(kc, kct, infior, i++);
}

/*
 * A first uniq marks a potential start of a region, a 2nd or later uniq, within
 * scope - keylength minus keywidth distance - triggers a region insertion or
 * update. Also retrieve for this index the offset.
 */
        // TODO: if all nextNts are the same increase keylength.

// reverse seq could be required for mapping - maybe 
// N-stretches, later also skips, future: SNVs, splice sites?
static int
ext_uq_bnd(kct_t *C kc, Hdr *C h)
{

    uint64_t p, dna = 0ul, rc = 0ul;   // first seq after skip
    uint64_t *kct = NULL, *k;
    unsigned i, index = 0, ext = kc->readlength - KEY_WIDTH;
    C uint32_t b2end = (*kc->bdit).e;
    C uint32_t b2start = (*kc->bdit).s;
    uint32_t b2pos = b2start;
    int res;

    // build 1st key
    for (p = b2pos + h->s_s; p != b2pos + h->s_s + KEY_WIDTH; ++p) {
        _EVAL(get_nextnt(kc, p));
        dna = _seq_next(res, dna, rc);
    }
    start_dbg(kc, h, dna);

    kc->kct_scope[0] = NULL;

    for (b2pos += KEY_WIDTH; b2pos < b2end; ++b2pos) { // until next uniq region, stretch or contig
        uint64_t nt;
        /* process new key */
        _get_new_kct(kc, kct, dna, rc);
        ++h->total;
EPQ0(dbg >5, "[%u%c]:\t", b2pos, IS_UQ(kct) ? '*' : ' ');
print_dna(dna, dbg >5);

        _EVAL(get_nextnt(kc, p++));
        if (IS_UQ(kct)) {
//print_dna(dna);

            if ((*kct & B2POS_MASK) == 0)
                *kct ^= p + 1;

        } else if (ALL_SAME_NTS(kct)) {
            // TODO: extension..?
            EPQ(IS_DBG_K(kc, kct), "[%u] got %lu", p, nt);

        } else {
            expect(CAN_ADD_ITER(kct)) {
                *kct += AT_NEXT;
            } else {
                *kct |= DISTINCT; // too many to verify, but safe and probably true.
            }
        }
        dna = _seq_next(res, dna, rc);
	if (IS_UQ(kct)) {

	    k = kc->kct_scope[0];
            if (k == NULL || (b2start + index + KEY_WIDTH - 1) == b2pos) {
                EPQ(dbg > 3, "(dnstrm adjoining uq %u)", h->mapable);
                decr_excise_all(kc, k, kct, index);
                (*kc->bdit).e = b2start;

            } else if (IS_UQ(k)) {
                //EPQ(dbg > 3, "(2nd+ uq %u)", h->mapable);

                if (((*kc->bdit).e + index) == b2pos) // a 2nd uq
                    h->mapable -= b2pos - index;

                decr_excise_all(kc, k, kct, index);

            } else {
                EPQ(dbg > 3, "(1st uq %u)", h->mapable);

                (*kc->bdit).e = b2pos;
                _EVAL(mark_all(kc, index, b2pos));
            }
	    index = 0;
	} else if (index == ext + 1) {
	    i = 0;
	    if (IS_UQ(kc->kct_scope[i])) {

                uint32_t last_uq_pos = b2pos - index;
                if ((*kc->bdit).e == b2start) {
                    (*kc->bdit).s = last_uq_pos - KEY_WIDTH + 1;
                    h->mapable += (*kc->bdit).s + 1 - b2start;
                    EPQ(dbg > 3, "(dnstrm adjoining, %u)", h->mapable);

                } else if (((*kc->bdit).e + index) == b2pos) {
                    EPQ(dbg > 3, "(postponed %u)", h->mapable);

                } else { // not downstream adjoining.

                    h->bnd.insert(kc->bdit, *kc->bdit); // insert a copy of the current
                    // start is earlier to enable building of the 1st key.
                    (*kc->bdit).s = last_uq_pos - KEY_WIDTH + 1;
                    h->mapable += last_uq_pos + 1;
                    EPQ(dbg > 3, "(kept %u)", h->mapable);

                }
		++i;
            }
	    _EVAL(mark_all(kc, index - i, b2pos));
	    index = 0;
	}
        ASSERT(index <= ext, return -EFAULT);
        kc->kct_scope[index++] = kct;
    }
    EPQ(IS_DBG_K(kc, kct), "last key was dbgtsoffs");
    if (b2pos != b2end) {
        EPR("b2pos != b2end: %u != %u (happens for Y on hg19", b2pos, b2end);
        show_mantras(kc, h);
        ASSERT(strncmp(kc->id + h->part[0], "Y", strlen("Y")) == 0, return -EFAULT);
    }
    k = kc->kct_scope[0];
    ASSERT(k != NULL, res = -EFAULT; goto err);
    if (IS_UQ(k)) {
        EPQ (dbg > 3, "Post loop adjoining boundary [%u, %u]", b2pos, h->mapable);
        if ((*kc->bdit).s == b2start){
            h->mapable = b2end - b2start - KEY_WIDTH;
EPQ(dbg > 3, "X:[%u, %u]:", b2pos, h->mapable);
            kc->bdit = h->bnd.erase(kc->bdit);
        } else if (((*kc->bdit).e + index) == b2pos) {
EPQ(dbg > 3, "E:[%u, %d]:", b2pos, h->mapable);
            h->mapable -= b2pos - index;
            h->mapable += b2pos + 1; // no extension beyond end
            kc->bdit++;
        } else {
            h->mapable += b2pos;

            kc->bdit++;
EPQ(dbg > 3, "last was uniq, %u", h->mapable);
            //show_mantras(kc, h);
        }
        decr_excise_all(kc, k, NULL, index);

    } else {
	_EVAL(mark_all(kc, index, b2pos));
        EPQ(dbg > 3, "(post_loop non-uq)");
//DESCRIBE_KEY(kc, k, 'k');
        (*kc->bdit).e = b2pos;
//        h->mapable += (*kc->bdit).s; //XXX
        ++kc->bdit;
        if (kc->bdit != h->bnd.end() && (*kc->bdit).e <= b2pos) {
            EPQ(dbg > 3, "FIXME: still needed?");
            //this seems necessary in case f N-stretch at end, why?
            //kc->bdit = h->bnd.erase(kc->bdit);
        }
    }
//show_mantras(kc, h);
    res = 0;
err:
    if (res < -1) {
        EPR0("following pos %u is: ", b2pos);
        for(i = 0, --p; i != KEY_WIDTH; ++i, ++p)
             EPR0("%c", b6(((kc->s[p>>2] >> ((p&3) << 1)) & 3) << 1));
        EPR0("\n[b2pos:%u, koffs:%lu], \t", b2pos, K_OFFS(kc, kct));
        print_2dna(dna, rc);
        show_mantras(kc, h);
    }
    return res;
}

static int
ext_uq_hdr(kct_t* kc, Hdr* h)
{
    //FIXME: rather than % of all display % mapable of remaining (mantras)
    kc->bdit = h->bnd.begin();
    h->mapable = h->total = 0ul;
    IFOUT(kc->bdit == h->bnd.end(), "Already fully mapable: %s", kc->id + h->part[0]);

    EPQ(dbg > 3, "Processing %s", kc->id + h->part[0]);

    do {
        int ret = ext_uq_bnd(kc, h);
        if (ret < 0) return ret;
    } while (kc->bdit != h->bnd.end());

    IFOUT(h->bnd.begin() == h->bnd.end(), "Became fully mapable: %s", kc->id + h->part[0]);
#ifdef DEBUG
    if (dbg > 4)
        show_mantras(kc, h);
#endif
    EPQ(dbg > 2, "Contig %s: %u/%u => %.2f%% mapable",
            kc->id + h->part[0], h->mapable, h->total,
            h->end_pos ? 100.0f * h->mapable / h->total : nanf("NAN"));
    ASSERT(h->mapable <= h->total, show_mantras(kc, h); return -EFAULT);
out:
    return 0;
}

static int
ext_uq_iter(kct_t* kc)
{
    kc->uqct = kc->pending = 0u;
    uint64_t mapable = 0ul, totNts = 0ul; // FIXME: put in kc and move to key_init
    for (std::list<Hdr*>::iterator h = kc->h.begin(); h != kc->h.end(); ++h) {
        // over ref headers
        int ret = ext_uq_hdr(kc, *h);
        if (ret < 0) return ret;
        mapable += (*h)->mapable;
        totNts += (*h)->total;
    }
    EPQ(dbg > 0, "extended %u unique ranges in iteration %u\n"
            "\t%lu/%lu => %.2f%% mapable. (%u pending)", kc->uqct, ++kc->iter,
            mapable, totNts, totNts ? 100.0f * mapable / totNts : nanf("NAN"), kc->pending);

    // It can occur that mapable != 0 while uqct == 0. This occurs when only adjoining uniques
    // are added in the iteration, but none are excised. In that case we stop looping nonetheless
    // because there cannot be any added new uniques in the subsequent iteration.
    return kc->uqct | kc->pending;
}

static int
extd_uniqbnd(kct_t* kc, struct gzfh_t* fhout)
{
    size_t ext = kc->readlength - KEY_WIDTH;
    int res = -ENOMEM;
    kc->iter = 0;

    kc->kct_scope = (uint64_t**)malloc((ext+1) * sizeof(uint64_t*));
    ASSERT(kc->kct_scope != NULL, res = -EFAULT; goto err);

    do { // until no no more new uniques
        res = ext_uq_iter(kc);
    } while (res > 0);
    _buf_free(kc->kct_scope);
    if (res == 0) {
        _ACTION(save_boundaries(fhout, kc), "writing unique boundaries file");
        _ACTION(save_kc(fhout + 3, kc), "writing unique keycounts file");
    }
err:
    return res;
}

int
fa_index(struct gzfh_t* fh, unsigned readlength)
{
    int len, res = -ENOMEM;
    char file[1024];
    kct_t kc = {0};
    kc.readlength = readlength;
    unsigned mode;

    C char* ext[7] = {".fa",  ".2b",".nn",".kc",".bd",  ".ub", ".uq"};

    if (fh[0].name) {
        len = strlen(fh[0].name);
    } else {
        ASSERT(fh[2].name != NULL, return -EFAULT);
        fh[0].name = file;

        len = strlen(fh[2].name);
        strncpy(fh[0].name, fh[2].name, ++len);
        _ACTION0(reopen(fh, ext[0], ext[5]), "")
    }
    ASSERT(len < 256, return -EFAULT, "filename too long: %s", fh[0].name);
    if (fh[1].name == NULL)
        fh[1].name = &file[len];

    strncpy(fh[1].name, fh[0].name, len);

    if (fh[3].name == NULL)
        fh[3].name = &file[len*2];

    strncpy(fh[3].name, fh[0].name, len);

    kc.ndxkct = _buf_init_arr_err(kc.ndxkct, KEYNT_BUFSZ_SHFT, return -ENOMEM);
    // first check whether unique boundary is ok.
    if (fh[0].fp) {
        mode = 2;
    } else {
         _ACTION0(reopen(fh, ext[5], ext[4]), "trying to open %s instead", ext[4])
        mode = fh[0].fp != NULL;
    }
    if (mode) {
        _ACTION(load_boundaries(fh, &kc), "loading boundary file %s", fh[0].name)
    }

    // keycount file did not exist, check twobit file.
    _ACTION0(reopen(fh, ext[mode == 2 ? 5 : 4], ext[1]),
                "%s does%s exist", ext[1], fh[0].fp ? "" : " not")
    _ACTION0(reopen(fh + 1, ext[5], ext[2]), "%s does%s exist", ext[2], fh[1].fp ? "" : " not")
    _ACTION0(reopen(fh + 3, ext[5], ext[3]), "%s does%s exist", ext[3], fh[3].fp ? "" : " not")

    if (mode && fh[0].fp && fh[1].fp && fh[3].fp) {
         _ACTION(load_seqb2(fh, &kc), "loading twobit sequence file")
         _ACTION(load_kc(fh + 3, &kc), "loading keycounts file")
    } else {
        bool found = false;
        for (int i=0; i != 4; ++i) {
            if (i == 2) continue;
            if (fh[i].fp) {
                found = true;
                EPR("%s file present, refusing to overwrite.", fh[i].name);
                rclose(fh +i);
            }
        }
        if (mode) {
            EPR("%s file present, refusing to overwrite.", ext[3+mode]);
            found = true;
        }
        if (found) goto err;
        EPR("starting from scratch.");
        _ACTION(fa_read(fh, &kc), "reading fasta")
    }
    if (mode < 2) {
        _ACTION(reopen(fh, ext[1], ext[5]), "")
        _ACTION(reopen(fh + 3, ext[3], ext[6]), "")
        _ACTION(extd_uniqbnd(&kc, fh), "extending unique boundaries")
    }
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free_kc(&kc);
    return res;
}


