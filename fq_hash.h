/*
 *    Description:
 *
 *        Version:  1.0
 *        Created:  16-06-14 21:17:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Roel KLuin,
 */

#ifndef FQ_ARR_H
#define FQ_ARR_H

#include <stdint.h>
#include <assert.h>
#include <stdio.h>

#include <iostream>
#include <algorithm>
#include <string.h>
#include "b6.h"
#include "util.h"

// maximum length of the read name
#define FQ_MAX_NAME_ETC    (1u << 14)

// KEY_LENGTH must be odd - or 2nd bit of central Nt is not always flipped in its
// complement - the alternative is the twisted halfdev conversion, but this is cheaper

#define KEY_LENGTH 15                       // <= this many Nts are used as key
#define KEY_CUTOFF 4                        // cut off, to get random keys
#define KEY_WIDTH (KEY_LENGTH + KEY_CUTOFF) // entire key maximized on (after conversion)

// KEYNT: 2bits, therefore multiplied by 2.
#define KEYNT_BUFSZ (1u << (KEY_LENGTH * 2 - 1))
#define KEYNT_TRUNC_MASK (KEYNT_BUFSZ - 1u)

#define KEYNT_STRAND (1ul << KEY_WIDTH)
#define KEYNT_AC (1ul << (KEY_WIDTH - 1))

#define KEYNT_TOP ((KEY_WIDTH - 1) * 2)
#define KEYNT_MASK ((1ul << KEY_WIDTH * 2) - 1ul)
#define HALF_KEYNT_MASK (KEYNT_STRAND - 1ul)

#define SEQ_OFFSET_BYTES 2
#define BUF_OFFSET_BYTES 4
#define SEQ_OFFSET_MAX ((1u << (SEQ_OFFSET_BYTES << 3)) - 1u)

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define if_ever(err)    if (unlikely(err))
#define expect(T)       if (likely(T))

// Exponentiated phred for summation R: as.integer(round(10^((0:63)/10)))
static unsigned phredtoe10[] = {
    1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 10, 13, 16, 20,
    25, 32, 40, 50, 63, 79, 100, 126, 158, 200, 251, 316, 398, 501,
    631, 794, 1000, 1259, 1585, 1995, 2512, 3162, 3981, 5012, 6310, 7943, 10000, 12589,
    15849, 19953, 25119, 31623, 39811, 50119, 63096, 79433, 100000, 125893, 158489, 199526,
    251189, 316228, 398107, 501187, 630957, 794328, 1000000, 1258925, 1584893, 1995262
};

template<typename T>
static inline void memset_T(void* dest, T value, uint64_t size)
{
  for(uint64_t i = 0; i != size; i += sizeof(T))
    memcpy(((char*)dest) + i, &value, sizeof(T));
}

using namespace std;

struct fq_arr
{
    fq_arr(unsigned rl, unsigned po) : is_err(0), l(0u), m(1u << 23),
        fq_ent_max((rl << 1) + FQ_MAX_NAME_ETC), nr(0u), readlength(rl),
        phred_offset(po), tidi(0u), tidm(32u)
    {
        s = (uint8_t*)malloc(m);
        *s = '\0';
        key_ct = 0u;
        tid = (unsigned*)malloc(tidm*sizeof(unsigned));
        unsigned looklen = KEYNT_BUFSZ;
        look = (uint32_t*)malloc(sizeof(uint32_t)*looklen);
        assert(look != NULL);
        memset_T<uint32_t>(look, 1u, sizeof(uint32_t)*looklen);
        kroundup32(fq_ent_max);
    }

    /**
     * Store or update the entry with the given sequence. The key is the gray maximum of
     * KEY_WIDTH Nts in the sequence, foreach the strand to be considered decided by the
     * central Nt's 2nd bit, see maximize_key(). As value is stored a 32 bit consisting
     * of an offset in buffer `s'.
     *
     * At the offset, buffer `s' holds two bytes, for the offset of the max in the sequence,
     * and a `complement state' in the highest bit of the fourth byte. Another four bytes point
     * to the offset in buffer `s' with stats for the former read that had this same deviant,
     * or `1' if there were none.
     *
     * Finally the sequence and phred quality are encoded in the buffer, ended by 0xff, followed
     * by its read name - '\0' terminated.
     *
     * For each next read with the same deviant key, the values' offset in buffer `s' is
     * updated to point to its last entry, until `1' which instead shows it was the last entry.
     */
    void put(const uint8_t* seq, const uint8_t* qual, const uint8_t* name)
    {
        ++nr;
        realloc_s();
        pre_init_key();
        if_ever (not init_key(seq)) return; // XXX: too short sequences are skipped.

        // get key - insensitive to strand. KEY_WIDTH Nts already read in init_key()
        uint8_t *b = add_max_key(seq + KEY_WIDTH);
        while ((*b++ = *name) != '\0') ++name;
        encode_seqphred(qual, seq, b);
        l = b - s;
    }
    /**
     * The key is the sequence or revcmp dependent on the 2nd bit of the central Nt. As every
     * 2nd bit of each twobit, it is inverted for the reverse complement. Since this is the
     * central Nt - odd KEY_LENGTH - the position also is the same for the reverse complement.
     *
     * For all windows of KEY_LENGTH plus KEY_CUTOFF, the key value is calculated. The maximum
     * thereof is kept as well as its offset.
     *
     */
    uint8_t* add_max_key(const uint8_t* sq)
    {
        uint64_t lmax, maxv, c;
        unsigned i;

        // seq or revcmp according to 2nd bit of central Nt
        c = (dna & KEYNT_STRAND) ? rev : dna;
        lmax = c;
        maxv = c ^ (c >> 2); // maximize on gray Nt sequence
        i = rot;

        while((c = *sq++) != '\0') { // search for maximum from left
            ++rot;
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            dna = ((dna << 2) & KEYNT_MASK & 0xffffffffffffffff) | c;
            rev = ((c ^ 2) << KEYNT_TOP) | (rev >> 2);
            c = (dna & KEYNT_STRAND) ? rev : dna;
            // with gray Nt code, we may select a more specific subsequence.
            uint64_t Ntgc = c ^ (c >> 2);

            if ((Ntgc > maxv) || unlikely(Ntgc == maxv && (dna & KEYNT_STRAND))) {
                lmax = c;
                maxv = Ntgc;
                i = rot;
            }
        }
        lmax = ((lmax >> 1) & ~HALF_KEYNT_MASK) | (lmax & HALF_KEYNT_MASK); // excise out 2nd cNt bit
        lmax &= KEYNT_TRUNC_MASK; /* truncate to make keys more random */

        // FIXME: assertion fails for large no. input reads.
        assert (l < (1ul << 32)); // or we cannot store its offset here
        uint8_t *b = s + l;
        // we could decrease this size to e.g. uint16_t: (jmax = 2);
        for (unsigned j = 0; j != SEQ_OFFSET_BYTES; ++j, i >>= 8) // insert key offset
            *b++ = i & 0xff;

        uint32_t* r = look + lmax;
        i = *r;

        if (i == 1u) key_ct++; // first element
        *r = l; // exchange by this offset

        for (unsigned j = 0; j != BUF_OFFSET_BYTES; ++j, i >>= 8)
            *b++ = i & 0xff; // move former buffer offset, or last element mark (1u) here
        return b;
    }
/*    void fa_put(const uint8_t* seq, unsigned* lmaxi)
    {
        realloc_s();
        if (cNt & 1) {// => last line was a new tid, so need to refill maxes
            if (cNt == 1) {i = 0; roti = -20; }
            cNt &= ~3;                           // unset `need-more-char bits'
            if (not init_key(seq)) { cNt |= 3; return; }
            key = encode_twisted(seq, &i);
            if (i == 0) { cNt |= 3; return; }    // need more chars
            seq += roti % (rl - 20);             // we already got the max for these
            update_offsets(key, i);
        }
        if (seq[0] != '>') {
            uint8_t *b = s + l;
            while((c = *sq++) != '\0') {
                if (++rot == rl - 20) rot = 0;

                c = b6(c);
                c = -isb6(c) & (c >> 1);
                dna = ((dna << 2) & KEYNT_MASK & 0xffffffffffffffff) | c;
                rev = ((c ^ 2) << KEYNT_TOP) | (rev >> 2);
                unsigned t = dna & KEYNT_STRAND;
                c = t ? rev : dna;
                // with gray Nt code, we may select a more specific subsequence.
                uint64_t Ntgc = c ^ (c >> 2);

                // rather than the max, we should search for the most distinct substring.
                if ((Ntgc > rot[*lmaxi]) || unlikely(Ntgc == rot[*lmaxi] && t)) {
                    *lmaxi = rot;
                    i = left_maxi | (-!(cNt ^ t) & 0x80000000);
                    key = c & 0xffffffff;
                    if_ever(insert_key(key, i) == false) // insert into database
                        return;

                } else if (roti == left_maxi) {
                    // the max has left, we have to search for the next max
                    // for each we also need to redetermine the twisted state
                    left_maxi = roti;
                    i = roti == 0 ? rl - 21 : roti - 1;
                    // if the second bit is the same for the last rotation
                    // the flip state was the same.
                    t = (t != 0);
                    t ^= (((rot[roti] >> 2) ^ rot[i]) & 0x2) != 0;
                    do {
                        // FIXME: calculate t:
                        if ((rot[i] > c) || unlikely(rot[i] == c && t)) {
                            left_maxi = i;
                        }
                        t ^= rot[i];
                        if (i + 1 == rl) {
                            t = rot[i] & 0x20000;
                            i = 0;
                    } while (i != roti);
                    c = rot[roti]
                }
            // if a max is leaving, search for the new max
            // insert every new char in buffer
            // for each new max, insert its key in the database
        } else if (seq[0] == '>') {
            ins_fa_hdr(seq);
            cNt = 1; // refill maxes and reset rotations
        }


        b = s + l;

        unsigned twisted, i = 0;
        // get twisted: value that is the same for seq and revcmp
        twisted = encode_twisted(seq, &i); // remove low bit, separate deviant
        if_ever (i == 0) return;           // XXX: means sequence too short. skipped.

        for (unsigned j = 0; j != 4; ++j, i >>= 8) // insert offset
            *++b = i & 0xff;

        // NOTE: the central Nt is not part of the key. A mismatch can occur at that Nt.
        khiter_t k = kh_get(UQCT, H, twisted);
        if (k == kh_end(H)) { // new deviant
            k = kh_put(UQCT, H, twisted, &is_err);
            if_ever ((is_err = (is_err < 0))) return;
            kh_val(H, k) = 1ul << RK_FQ_KEYCOUNT_OFFSET; // keycount
            i = 1; // marks last element

        } else {
            kh_val(H, k) += 1ul << RK_FQ_KEYCOUNT_OFFSET;

            i = kh_val(H, k) & 0xffffffff; // former value offset
            kh_val(H, k) ^= i;
        }
        kh_val(H, k) |= l;

        for (unsigned j = 0; j != 4; ++j, i >>= 8) // last insert offset
            *++b = i & 0xff;

        b = encode_seqphred(qual, seq, b);
        while ((*++b = *name) != '\0') ++name;
        l = ++b - s;
     }*/

    /**
     * print the fastq entries. 
     */
    void dump()
    { /* must always be called to clean up */
        if (not is_err) {
            uint8_t* d = (uint8_t*)malloc((readlength<<1)+3);
            char seqrc[(KEY_LENGTH+1)<<1];
            cout << hex;

            uint32_t j = KEYNT_BUFSZ - 1u;
            do {
                uint32_t i = look[j];
                if (i != 1u) { // has keys
                    decode_key(seqrc, j);
                    do {
                        uint8_t* b = s + i + SEQ_OFFSET_BYTES + BUF_OFFSET_BYTES - 1;
                        i = *b; i <<= 8;
                        i |= *--b; i <<= 8;
                        i |= *--b; i <<= 8;
                        i |= *--b; // next seq offset for same key, or 1u.
                        // similarly key offset could be retrieved, but KEY_WIDTH
                        // needs to be added to get to the real offset.
                        b += BUF_OFFSET_BYTES;
                        fprintf(stdout, "%s 0x%x\t%s\n", (const char*)s, j, seqrc);
                        decode_seqphred(b + strlen((const char*)s) + 1, d);

                        fprintf(stdout, "%s\n+\n%s\n", d, d + readlength + 1);
                    } while (i != 1u);
                }
            } while (--j != 0u);
            free(d);
        }
        free(tid);
	free(s);
        free(look);
    }
    void pre_init_key() { dna = rev = 0ul; rot = KEY_WIDTH; }
    bool init_key(const uint8_t* sq)
    {
        do {
            uint64_t c;
            if ((c = *sq++) == '\0') return false; // need more sequence 
            // zero [min] for A or non-Nt. Maybe reads with N's should be processed at the end.
            c = b6(c);
            c = -isb6(c) & (c >> 1);
            dna = (dna << 2) | c;
            rev = (c << KEYNT_TOP) | (rev >> 2);
        } while (--rot != 0u);
        rev ^= KEYNT_MASK & 0xaaaaaaaaaaaaaaaa; // make reverse complement
        return true;
    }

    int is_err;

private:
    void realloc_s()
    {   assert (l != 1u);
        assert((uint64_t)l + fq_ent_max < (1ul << 32));
        if (l + fq_ent_max >= m) { // grow buffer if insufficient space for another read
            m <<= 1;
            s = (uint8_t*)realloc(s, m);
        }
    }
    void ins_fa_hdr(const uint8_t* seq)
    {
        uint8_t *b = s + l;
        // insert position of next reference string.
        if (tidi++ == tidm) {
            tidm <<= 1;
            tid = (unsigned*)realloc(tid, tidm*sizeof(unsigned));
        }
        *tid++ = l;

        // insert reference string, skip '>'
        while ((*b++ = *++seq) != '\0') {}
        l = b - s;
    }

    void decode_key(char* seq, unsigned dna)
    {
        char* revcmp = seq + (KEY_LENGTH << 1) + 1;
        *revcmp = '\0';
        //put back cNT bit in seqrc (zero)
        dna = ((dna << 1) & ~HALF_KEYNT_MASK) | (dna & HALF_KEYNT_MASK);

        char c;

        for (unsigned i = 0; i != KEY_LENGTH; ++i, ++seq, dna >>= 2) {
            c = (dna & 0x3) << 1;
            *seq = b6(c ^ 4);
            *--revcmp = b6(c);
        }
        *seq = '|';
        assert(revcmp == seq + 1);
    }

    /**
     * put sequence and phred qualities in one char foreach Nt. This includes the key
     * sequence section currently, but this may change in the future - it is not needed.
     */
    void encode_seqphred(const uint8_t* qual, const uint8_t* seq, uint8_t *b)
    {
        unsigned i = 0;
        while ((*b = *qual++) != '\0') {
            *b -= phred_offset;
            assert(*b <= 50);
            i += phredtoe10[*b];
            unsigned c = b6(*seq++);

            // In case of an 'N' phred should always be less than 3
            // max phred for base (G) is 60 using this scheme
            if (isb6(c)) *b = (*b | (c << 5)) + 3;
            ++b;
        }
        *b = 0xff;
        l = ++b - s;
    }
    /**
     * decode the sequence and phreds.
     */
    void decode_seqphred(uint8_t* b, uint8_t* d)
    {
        uint8_t *r = d + readlength + 1;

        for(; *b != 0xff; ++b, ++d, ++r) {
            unsigned c = *b;
            if (c >= 3) {
                c -= 3;
                *d = b6((c >> 5) & 0x6);
                *r = (c & 0x3f) + phred_offset;
            } else {
                *d = 'N';
                *r = c + phred_offset;
            }
        }
        *d = *r = '\0';
    }
    uint32_t l, m, fq_ent_max, nr, readlength;
    unsigned phred_offset, key_ct, rot, tidi, tidm;
    uint8_t *s;
    unsigned* tid;
    uint32_t* look;
    uint64_t dna, rev;
};


#endif //FQ_ARR_H

