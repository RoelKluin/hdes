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
#include "fa.h"

/**
 * Initialize buffer and fill it with 1's. Safe because the first entry starts at 0.
 */
static inline int
init_fq(seqb2_t *seq)
{

    seq->s = _buf_init(seq->s, INIT_BUFSIZEBIT);
    *seq->s = '\0';

    size_t i, v = 1, l = sizeof(*seq->lookup) * KEYNT_BUFSZ;
    seq->lookup = (uint32_t*)malloc(l);
    if (seq->lookup == NULL) {
        _buf_free(seq->s);
        return -ENOMEM;
    }
    for (i = 0; i != l; i += sizeof(*seq->lookup))
        memcpy(((char*)seq->lookup) + i, &v, sizeof(*seq->lookup));
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
get_tid_and_pos(kct_t* kc, uint64_t *pos, unsigned bufi)
{
    // FIXME: no looping here, store ref to header in boundary?.
    std::list<Hdr*>::iterator hdr;
    std::list<uint32_t>::iterator bd;
    for (hdr = kc->h.begin(); hdr != kc->h.end(); ++hdr) {
        bd = (*hdr)->bnd.end();
        --bd;
        ASSERT(*bd < kc->bd_l, return -EFAULT);
        //position beyond end of last boundary of this contig must be on a later one.
        if (((*hdr)->s_s + kc->bd[*bd].s + kc->bd[*bd].l) < *pos)
            continue;
        EPQ(dbg > 16, "%s", kc->id + (*hdr)->part[0]);
        while ((((*hdr)->s_s + kc->bd[*bd].s + kc->bd[*bd].l) > *pos) &&
                bd != (*hdr)->bnd.begin()) {
            EPQ(dbg > 16, "%lu > %lu?", ((*hdr)->s_s + kc->bd[*bd].s + kc->bd[*bd].l), pos);
            --bd;
        }
        ASSERT (((*hdr)->s_s + kc->bd[*bd].s + kc->bd[*bd].l) <= *pos, return -EFAULT);
        break;
    }
    ASSERT (hdr != kc->h.end(), return -EFAULT);
    EPQ(dbg > 6, "%ld <= %u + %u", *pos - (*hdr)->s_s, kc->bd[*bd].corr, bufi);

    ASSERT(*pos + kc->bd[*bd].corr > (*hdr)->s_s + bufi, return -EFAULT/*,
            "\n%lx + %x <= %lx + %x +s", *pos, kc->bd[*bd].corr, (*hdr)->s_s, bufi, KEY_WIDTH*/)
    *pos += kc->bd[*bd].corr - (*hdr)->s_s - bufi; // should be one-based.

    return (*hdr)->part[0];
}


static int
fq_read(kct_t* kc, seqb2_t *seq)
{
    void* g;
    int (*gc) (void*);
    int (*ungc) (int, void*);
    struct gzfh_t* fhin = seq->fh;
//dbg = 7;
    uint64_t l = seq->s_l, m = seq->s_m;
    register uint8_t *s = seq->s + l;
    const unsigned phred_offset = seq->phred_offset;
    unsigned fq_ent_max = SEQ_MAX_NAME_ETC + seq->readlength + 1;
    register int c;
    uint64_t ndx, dna = 0ul, rc = 0ul, b = 0ul, wx = 0ul;
    uint64_t* buf = (uint64_t*)malloc((kc->ext + 1) * sizeof(uint64_t));
    unsigned* bufi = (unsigned*)malloc((kc->ext + 1) * sizeof(unsigned));
    Hdr* lh = kc->h.back();
    const uint64_t end_pos = lh->s_s + lh->end_pos;

    set_readfunc(fhin, &g, &gc, &ungc);

    while ((c = gc(g)) != '@') /* skip to first header */
        if (c == -1 || c == '>') goto out;
    do {
        register unsigned i = 0;
        _buf_grow0(seq->s, fq_ent_max); // at least enough space for one read

        s = seq->s + seq->s_l;
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

        char* seqstart = (char*)s;
        i = 0;
        while ((c = gc(g)) != '\n') {
            ASSERT(c != -1, c = -EFAULT; goto out);
            *s = b ^= b;
            switch(c) {
case 'C': case 'c': b = 2;
case 'G': case 'g': b ^= 1;
case 'U': case 'u':
case 'T': case 't': b ^= 2;
case 'A': case 'a': *s = 0x3e;
default:            dna = _seq_next(b, dna, rc);
                    *s++ |= (b << 6) | 1; // for seqphred storage
                    if (++i < KEY_WIDTH) // first only complete key
                        continue;
		    ASSERT(i <= kc->ext + KEY_WIDTH, c = -EFAULT; goto out);
		    // wx only has strand at offset KEY_WIDTH.
                    _get_ndx(ndx, wx, dna, rc);
                    uint32_t k = kc->ndxkct[ndx];
		    // put not recognized and multimapper keys to end - unused.
                    if (k >= kc->kct_l ||
                            (kc->kct[k + 1] >> BIG_SHFT) > 1ul) {
                        buf[i - KEY_WIDTH] = ~0ul; // FIXME: could write ndx here.
                        bufi[i - KEY_WIDTH] = i ^ wx;
                        continue;
                    }
                    ASSERT((kc->kct[k] & B2POS_MASK) < end_pos,
                        c = -EFAULT; goto out, "0x%lx\t%d", ndx, print_ndx(ndx));
                    ASSERT((int)(kc->kct[k] & B2POS_MASK) >= 0,
                            c = -EFAULT; goto out, "ndx:0x%lx\n[0]:0x%lx\n[1]:0x%lx\npos:%d", ndx,
                            kc->kct[k],
                                kc->kct[k + 1],
                            (int)(kc->kct[k] & B2POS_MASK) & print_ndx(ndx));
                    ASSERT(k != -2u, c = -EFAULT; goto out)
                    if (dbg > 6) {
                        uint64_t pos = kc->kct[k] & B2POS_MASK;
                        c = get_tid_and_pos(kc, &pos, i);
			if (c < 0) continue;
                        ASSERT(c >= 0, goto out, "%lu", kc->kct[k] & B2POS_MASK);
                        EPR0("%s:%lu\t%lu\t0x%lx\t", kc->id + c, pos,
                                kc->kct[k] >> INFIOR_SHFT, ndx);
                        print_ndx(ndx);
                    }
		    if (i == KEY_WIDTH) {
	                    buf[0] = k;
			    bufi[0] = i ^ wx;
		    } else if (buf[0] == ~0ul){
			buf[i - KEY_WIDTH] = buf[0];
			bufi[i - KEY_WIDTH] = bufi[0];
			buf[0] = k;
			bufi[0] = i ^ wx;
		    } else {
			// test infior - in high bits.
                        uint32_t t = buf[0];
                        // TODO: early verify and process unique count if correct.
                        ASSERT(k != -2u, c = -EFAULT; goto out);
                        if ( kc->kct[k] > kc->kct[t]) {
                            buf[i - KEY_WIDTH] = kc->ndxkct[ndx];
		            bufi[i - KEY_WIDTH] = i ^ wx;
			} else {
			    // lowest inferiority
                            buf[i - KEY_WIDTH] = buf[0];
                            bufi[i - KEY_WIDTH] = bufi[0];
                            buf[0] = k;
                            bufi[0] = i ^ wx;
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

        ASSERT((kc->kct[ndx] & B2POS_MASK) < end_pos,
                c = -EFAULT; goto out, "0x%lx", ndx);
        //ASSERT((int)(kc->kct[ndx] & B2POS_MASK) >= 0,
        //        c = -EFAULT; goto out, "0x%lx*\n0x%lx\n%d", ndx, 
        //                    kc->kct[ndx],
        //                    (int)(kc->kct[ndx] & B2POS_MASK));

        if (((uint32_t)ndx < kc->kct_l) &&
                (kc->kct[ndx + 1] >> BIG_SHFT) == 1ul) {
            mq = 37;
            flag = (wx ^ (kc->kct[ndx] >> BIG_SHFT)) & 1;
            if ((wx & 1)) { //XXX
                while ((c = gc(g)) != -1 && c != '@') {}
                continue;
            }
	    flag = !(kc->kct[ndx] & STRAND_BIT) ^ !(bufi[0] & KEYNT_STRAND);
            //ASSERT(flag == 0, return -EFAULT, "FIXME: ASSERTION only valid in testset");
            flag = 0x42 | (flag << 4);



            uint64_t pos = kc->kct[ndx] & B2POS_MASK;

            EPQ(dbg > 6, "pos:%lu", pos);
            c = get_tid_and_pos(kc, &pos, bufi[0] & ~KEYNT_STRAND);
	    if (c < 0) {	
		while ((c = gc(g)) != '@' && c != -1) {}
		continue;
	    }
            ASSERT(c >= 0, goto out, "%lu", kc->kct[ndx] & B2POS_MASK);
            ASSERT((int)pos > 0, c = -EFAULT; goto out,
                    "pos:%lu, kc->kct[0x%lx,+1], %lx, %lx", pos, ndx, kc->kct[ndx], kc->kct[ndx+1])

            const char* tid = kc->id + c;

            OPR0("%s\t%u\t%s\t%lu\t%u\t%uM\t%s\t%u\t%d\t",
           (char*)h,flag,tid, pos,mq,seqlen,mtd,mps,tln);
            if (flag & 0x10) { //reverse complemented
                for (i = 0; i != seqlen; ++i) {
                    c = *--s;
                    *s = gc(g);
                    ASSERT(*s != '\n' && *s != -1, c = -EFAULT; goto out, "truncated");
                    switch (c) {
                        case 0x3f: OPR0("T"); break;
                        case 0x7f: OPR0("G"); break;
                        case 0xbf: OPR0("A"); break;
                        case 0xff: OPR0("C"); break;
                        case 0x1: OPR0("N"); break;
                    }
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
    free(buf);
    free(bufi);
    return c;
}

int
map_fq_se(struct seqb2_t* seq)
{
    int res = -ENOMEM;

    // 1) open keyindex, infior and strand
    struct gzfh_t* fhio[3] = { seq->fh + 1, seq->fh + 2, seq->fh + 3};
    const char* ext[4] = {".kc",".2b",".bd",  ".uq"};
    char file[768];
    kct_t kc = {0};
    kc.ext = seq->readlength - KEY_WIDTH;
    ASSERT(fhio[1]->name != NULL, return -EFAULT);
    unsigned len = strlen(fhio[1]->name) + 1;
    ASSERT(strstr(fhio[1]->name, ext[1]), return -EFAULT);

    _ACTION( init_fq(seq), "intializing memory")

    for (int i=0; i != 3; ++i) {
        if (fhio[i]->name == NULL) {
            fhio[i]->name = &file[len*i];
            strncpy(fhio[i]->name, fhio[1]->name, len);
            _ACTION0(reopen(fhio[i], ext[1], ext[i]), 
                    "%s %s", ext[i], fhio[i]->fp ? "exists" : "does not exist")
        }
    }
    ASSERT(fhio[0]->fp && fhio[1]->fp && fhio[2]->fp, return -EFAULT,
            "need seqb2, keycount and unique boundary files");
    kc.ndxkct = _buf_init_arr_err(kc.ndxkct, KEYNT_BUFSZ_SHFT, return -ENOMEM);
    // 2) open seqb2 for verification of reads
    // 3) open original boundaries
    _ACTION(load_kc(fhio[0], &kc), "loading keycounts file")
    _ACTION(load_seqb2(fhio[1], &kc), "loading twobit sequence file")
    _ACTION(load_boundaries(fhio[2], &kc), "loading boundary file")

    _ACTION(reopen(fhio[0], ext[0], ext[3]), "")
    _ACTION(ammend_kc(fhio[0], &kc), "ammending keycounts from file")
    
    // 4) open fq for reading
    _ACTION(fq_read(&kc, seq), "mapping reads")
    //...
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free_kc(&kc);
    return res;
}
