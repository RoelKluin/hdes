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
#include <errno.h> // ENOMEM
#include <string.h> // memset()
#include "map.h"

/**
 * Initialize buffer and fill it with 1's. Safe because the first entry starts at 0.
 */
static inline int
init_fq(seqb2_t *sb2)
{

    sb2->s = _buf_init(sb2->s, INIT_BUFSIZEBIT);
    *sb2->s = '\0';

    size_t i, v = 1, l = sizeof(*sb2->lookup) * KEYNT_BUFSZ;
    sb2->lookup = (uint32_t*)malloc(l);
    if (sb2->lookup == NULL) {
        _buf_free(sb2->s);
        return -ENOMEM;
    }
    for (i = 0; i != l; i += sizeof(*sb2->lookup))
        memcpy(((char*)sb2->lookup) + i, &v, sizeof(*sb2->lookup));
    return 0;
}

/*
 * store seqphred and return corresponding 2bit Nt
 */
static inline unsigned
seqphred(uint8_t *s, int q)
{
    unsigned c = b6(*s);
    *s = q;

    // 1) seqphred includes redundant key Nts. We could save a few bytes per sequence.
    // 2) Maybe reads with N's should be processed at the end - in assembly.
    if (isb6(c)) {
        *s = (*s | (c << 5)) + 3;
        return c >>= 1;
    }
    return c ^ c; // zero [min] for N - but also for A.
}


static int
get_tid_and_pos(kct_t* kc, uint64_t *pos, C unsigned bufi)
{
    // FIXME: no looping here, store ref to header in boundary?.
    std::list<Mantra>::iterator bd;
    Hdr* h;

    for (h = kc->h; h != kc->h + kc->h_l; ++h) {
        bd = h->bnd->end();
        --bd;

        //position beyond end of last boundary of this contig must be on a later one.
        if ((h->s_s + (*bd).e + bufi) < *pos)
            continue;

        EPQ(dbg > 16, "%s", kc->id + h->part[0]);
        while (((h->s_s + (*bd).s + bufi) > *pos) && bd != h->bnd->begin()) {

            EPQ(dbg > 16, "%lu > %lu?", (h->s_s + (*bd).s), *pos);
            --bd;
        }
        ASSERT ((h->s_s + (*bd).s) <= *pos, return -EFAULT);
        break;
    }
    ASSERT (h != kc->h + kc->h_l, return -EFAULT);
    EPQ(dbg > 6, "%ld <= %u + %u", *pos - h->s_s, (*bd).corr, bufi);

    ASSERT(*pos + (*bd).corr > h->s_s + bufi, return -EFAULT,
            "%s\t%lu\t%lu", kc->id + h->part[0], *pos + (*bd).corr, h->s_s + bufi
            /*, "\n%lx + %x <= %lx + %x +s", *pos, (*bd).corr, h->s_s, bufi, KEY_WIDTH*/);
    *pos += (*bd).corr - h->s_s + KEY_WIDTH - bufi - 1;

    return h->part[0];
}

static void
print_hdr(kct_t *C kc, char C*C commandline)
{
    std::list<Hdr*>::iterator h;
    for (Hdr* h = kc->h; h != kc->h + kc->h_l; ++h) {
        OPR("@SQ\tSN:%s\tLN:%u", kc->id + h->part[0], h->end_pos + h->bnd->back().corr);
        EPR("@SQ\tSN:%s\tLN:%u", kc->id + h->part[0], h->end_pos + h->bnd->back().corr);
    }
    OPR("@PG\tID:" PROGRAM_NAME "\tPN:" PROGRAM_NAME "\tVN:" PROGRAM_VERSION "\tCL:%s",
            commandline);
}


static int
fq_read(kct_t* kc, seqb2_t *sb2)
{
    void* g;
    int (*gc) (void*);
    struct gzfh_t* fhin = sb2->fh;
//dbg = 7;
    uint64_t l = sb2->s_l;//, m = sb2->s_m;
    uint8_t *s = sb2->s + l;
    //const unsigned phred_offset = sb2->phred_offset;
    unsigned fq_ent_max = SEQ_MAX_NAME_ETC + sb2->readlength + 1;
    seq_t ndx;
    int c = kc->readlength - KEY_WIDTH + 1;
    uint64_t* buf = (uint64_t*)malloc(c * sizeof(uint64_t));
    unsigned* bufi = (unsigned*)malloc(c * sizeof(unsigned));
    keyseq_t seq = {0};
    Hdr* lh = kc->h + kc->h_l - 1;
    const uint64_t end_pos = lh->s_s + lh->end_pos;
    //struct mapstat_t ms = {0};

    set_readfunc(fhin, &g, &gc);

    while ((c = gc(g)) != '@') {/* skip to first header */
        ASSERT(c != -1 && c != '>', goto out);
    }
    do {
        unsigned i = 0;
        _buf_grow0(sb2->s, fq_ent_max); // at least enough space for one read

        s = sb2->s + sb2->s_l;
        uint8_t* h = s;

        while (!isspace(c = gc(g)) && (c >= 0)) *s++ = c; /* header */
        while (c != '\n' && (c >= 0)) c = gc(g);          /* comment - ignored. */
        if (c == -1) break;
        *s++ = '\0';
        if (strcmp((const char*)h, dbgrn)) {
            dbg &= ~8;
        } else {
            dbg |= 8;
        }
        EPQ(dbg > 6, "%s", (char*)h);

        //char* seqstart = (char*)s;
        i = 0;
        while ((c = gc(g)) != '\n') {
            ASSERT(c != -1, c = -EFAULT; goto out);
            *s = seq.t = 0;
            switch(c) {
case 'C': case 'c': seq.t = 2;
case 'G': case 'g': seq.t ^= 1;
case 'U': case 'u':
case 'T': case 't': seq.t ^= 2;
case 'A': case 'a': *s = 0x3e;
default:            seq_next(seq);
                    *s++ |= (seq.t << 6) | 1; // for seqphred storage
                    if (++i < KEY_WIDTH) // first only complete key
                        continue;
		    if (i > kc->readlength) { // none of keys mappable
			EPR("Read %s is with %u Nts longer than expected (%u) and skipped.",
				(char*)h, i, kc->readlength);
			buf[0] = ~0ul; //skip entire read
			break;
		    }
		    //ASSERT(i <= kc->readlength, c = -EFAULT; goto out);
		    // seq.t only has strand at offset KEY_WIDTH.
                    ndx = get_ndx(seq);
                    uint32_t k = kc->ndxkct[ndx];
 EPR("%u:%lx\t%x", i, (uint64_t)ndx, k);
		    // put not recognized and multimapper keys to end - unused.
                    if (k >= kc->kct_l || IS_UQ(kc->kct + k) == false) {
                        buf[i - KEY_WIDTH] = ~0ul; // FIXME: could write ndx here.
                        bufi[i - KEY_WIDTH] = i ^ seq.t;
                        continue;
                    }
		    if (b2pos_of(kc->kct[k]) >= end_pos) { // beyond chromosomes?
                        buf[i - KEY_WIDTH] = ~0ul; // FIXME: could write ndx here.
                        bufi[i - KEY_WIDTH] = i ^ seq.t;
                        continue;
                    }
                    if (dbg > 6) {
                        uint64_t pos = b2pos_of(kc->kct[k]);
                        c = get_tid_and_pos(kc, &pos, i);
			if (c < 0) {
                            buf[i - KEY_WIDTH] = ~0ul;
                            bufi[i - KEY_WIDTH] = i ^ seq.t;
                            continue;
                        }
                        EPR0("%s:%lu\t%lu\t0x%lx\t", kc->id + c, pos,
                                kc->kct[k] >> INFIOR_SHFT, (uint64_t)ndx);
                        print_ndx(ndx);
                    }
		    if (i == KEY_WIDTH) {
	                    buf[0] = k;
			    bufi[0] = i ^ seq.t;
		    } else if (buf[0] == ~0ul){
			buf[i - KEY_WIDTH] = buf[0];
			bufi[i - KEY_WIDTH] = bufi[0];
			buf[0] = k;
			bufi[0] = i ^ seq.t;
		    } else {
			// test infior - in high bits.
                        uint32_t t = buf[0];
                        // TODO: early verify and process unique count if correct.
                        if (kc->kct[k] > kc->kct[t]) {
                            buf[i - KEY_WIDTH] = kc->ndxkct[ndx];
		            bufi[i - KEY_WIDTH] = i ^ seq.t;
			} else {
			    // lowest inferiority
                            buf[i - KEY_WIDTH] = buf[0];
                            bufi[i - KEY_WIDTH] = bufi[0];
                            buf[0] = k;
                            bufi[0] = i ^ seq.t;
                        }
                    }
            }
        }
        if (buf[0] == ~0ul) {
            //EPR("FIXME: skip too short read or handle entirely non-matching");
            while ((c = gc(g)) != -1 && c != '@') {}
            continue;
        }
        ndx = buf[0];
        *s = '\0';
        while ((c = gc(g)) != '\n' && c != -1) {} // skip 2nd hdr line
        unsigned seqlen = i, tln = 0, mps = 0, mq = 0, flag = 44;
        const char* mtd = "*";
	if (b2pos_of(kc->kct[ndx]) >= end_pos) {
            while ((c = gc(g)) != -1 && c != '@') {}
            continue;
        }



        if (((uint32_t)ndx < kc->kct_l) && IS_UQ(kc->kct + ndx)) {
            mq = 37;
            flag = (seq.t ^ 1) & 1;
            if ((seq.t & 1)) { //XXX
                while ((c = gc(g)) != -1 && c != '@') {}
                continue;
            }
	    flag = !(kc->kct[ndx] & STRAND_BIT) ^ !(bufi[0] & KEYNT_STRAND);
            flag = (flag << 4) | 0x42;



            uint64_t pos = b2pos_of(kc->kct[ndx]);

            EPQ(dbg >> 6, "pos:%lu", pos);
            c = get_tid_and_pos(kc, &pos, bufi[0] & ~KEYNT_STRAND);
	    if (c < 0) {
		while ((c = gc(g)) != '@' && c != -1) {}
		continue;
	    }
            const char* tid = kc->id + c;

            OPR0("%s\t%u\t%s\t%lu\t%u\t%uM\t%s\t%u\t%d\t",
           (char*)h,flag,tid, pos,mq,seqlen,mtd,mps,tln);
            if (flag & 0x10) { //reverse complemented
                for (i = 0; i != seqlen; ++i) {
                    c = *--s;
                    switch (c) {
                        case 0x3f: OPR0("T"); break;
                        case 0x7f: OPR0("G"); break;
                        case 0xbf: OPR0("A"); break;
                        case 0xff: OPR0("C"); break;
                        case 0x1: OPR0("N"); break;
                    }
                    c = gc(g);
                    ASSERT(c != '\n' && c != -1, goto out, "truncated");
                    *s = c;
                }
                 OPR0("\t");
                for (i = 0; i != seqlen; ++i)
                    OPR0("%c", *s++);
                ASSERT(gc(g) == '\n', c = -EFAULT; goto out);
            } else {
                s -= seqlen;
                for (i = 0; i != seqlen; ++i) {
                    switch ((c = *s++)) {
                        case 0x3f: OPR0("A"); break;
                        case 0x7f: OPR0("C"); break;
                        case 0xbf: OPR0("T"); break;
                        case 0xff: OPR0("G"); break;
                        case 0x1: OPR0("N"); break;
                    }
                }
                OPR0("\t");
                while ((c = gc(g)) != '\n' && c != -1)
                    OPR0("%c", c);
            }
            OPR("\tMD:Z:%u\tNM:i:0", seqlen);
            s = h; // no storage needed
        } else {
            // XXX: skip for now.
            while ((c = gc(g)) != '\n' && c != -1) {}
            s = h; // no storage needed
        }
        c = gc(g);
    } while (c >= 0);
    c = 0;
out:
    QARN(c < 0, "error:%d", c);
    free(buf);
    free(bufi);
    return c;
}

int
map_fq_se(struct seqb2_t* sb2, char C*C cmdl)
{
    int res = -ENOMEM;

    // 1) open keyindex, infior and strand
    struct gzfh_t* fhio[3] = { sb2->fh + 1, sb2->fh + 2, sb2->fh + 3};
    const char* ext[4] = {".kc",".2b",".bd",  ".uq"};
    char file[768];
    kct_t kc = {0};
    kc.readlength = sb2->readlength;
    ASSERT(fhio[1]->name != NULL, return -EFAULT);
    unsigned len = strlen(fhio[1]->name) + 1;
    ASSERT(strstr(fhio[1]->name, ext[1]), return -EFAULT);
    // TODO: first read in several reads to verify 

    _ACTION(init_fq(sb2), "intializing memory");

    for (int i=0; i != 3; ++i) {
        if (fhio[i]->name == NULL) {
            fhio[i]->name = &file[len*i];
            strncpy(fhio[i]->name, fhio[1]->name, len);
            _ACTION0(reopen(fhio[i], ext[1], ext[i]), 
                    "%s %s", ext[i], fhio[i]->fp ? "exists" : "does not exist");
        }
    }
    ASSERT(fhio[0]->fp && fhio[1]->fp && fhio[2]->fp, return -EFAULT,
            "need seqb2, keycount and unique boundary files");
    kc.ndxkct = _buf_init_arr_err(kc.ndxkct, KEYNT_BUFSZ_SHFT, return -ENOMEM);
    // 2) open seqb2 for verification of reads
    // 3) open original boundaries
    _ACTION(load_kc(fhio[0], &kc), "loading keycounts file");
    _ACTION(load_seqb2(fhio[1], &kc), "loading twobit sequence file");
    _ACTION(load_boundaries(fhio[2], &kc), "loading boundary file");

    _ACTION(reopen(fhio[0], ext[0], ext[3]), "");
    _ACTION(ammend_kc(fhio[0], &kc), "ammending keycounts from file");
    
    // 4) print header
    print_hdr(&kc, cmdl);
    // 5) open fq for reading
    _ACTION(fq_read(&kc, sb2), "mapping reads");
    //...
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free(sb2->lookup);
    free(sb2->s);
    free_kc(&kc);
    return res;
}
