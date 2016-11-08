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
#include <ctype.h> // isspace() toupper()
#include <errno.h> // ENOMEM
#include <string.h> // memset()
#include "map.h"

/**
 * Initialize buffer and fill it with 1's. Safe because the first entry starts at 0.
 */
static inline int
init_fq(seqb2_t *sb2)
{

    sb2->s = buf_init(sb2->s, INIT_BUFSIZEBIT);
    *sb2->s = '\0';

    size_t i, v = 1, l = sizeof(*sb2->lookup) * KEYNT_BUFSZ;
    sb2->lookup = (uint32_t*)malloc(l);
    if (sb2->lookup == NULL) {
        buf_free(sb2->s);
        return -ENOMEM;
    }

    for (i = 0; i != l; i += sizeof(*sb2->lookup))
        memcpy(((char*)sb2->lookup) + i, &v, sizeof(*sb2->lookup));

    return 0;
}



static void
print_hdr_lines(Key_t *C kc, char C*C commandline)
{
    for (Hdr* h = kc->h; h != kc->h + kc->h_l; ++h) {
        OPR("@SQ\tSN:%s\tLN:%u", kc->id + h->ido, h->end + h->corr);
        EPR("@SQ\tSN:%s\tLN:%u", kc->id + h->ido, h->end + h->corr);
    }
    OPR("@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s", PROGRAM_NAME, PROGRAM_NAME,
            PROGRAM_VERSION, commandline);
}

static inline void
skip_plus_line(gzin_t* gz)
{
    int c;
    NB(GZTC(gz) == '\n');
    NB(GZTC(gz) == '+');
    while ((c = GZTC(gz)) != '\n')
        NB(c != -1);
}

static inline void
print_hdr(gzin_t* gz)
{
    int c;
    NB(GZTC(gz) == '@');

    while ((c = GZTC(gz)) != '\n') {
        NB(c != -1);
        fputc(c, stdout);
    }
    fputc('\t', stdout);
}

static int
match_revcmp(gzin_t* gz, struct map_t &map)
{
    // There can still be mismatches in upseq or post key.
    unsigned mismatches = 0;

    uint32_t p = (map.p >> 1) + map.iend - 1;
    EPR("p:%u", p);
    uint8_t* s = map.s + (p >> 2);
    char* upseq = map.upseq;
    unsigned ref_nt = *s >> ((p & 3) << 1);

    for (unsigned j = 0; j != map.readlength; ++j) {

        if (j >= map.iend) {
            int c = GZTC(gz);
            if (c == '\n' || c == -1)
                return -1;
            *upseq = toupper(c);
        }

        unsigned read_nt = 2;
        switch(*upseq++) {
            case 'C': read_nt = 0;
            case 'G': read_nt |= 1;
            case 'U':
            case 'T': read_nt ^= 2;
        }
        EPR("%c <-> %c", (char)b6(read_nt << 1), (char)b6(ref_nt << 1));

        if (read_nt != (ref_nt & 3) ) {

            // TODO: mism[j >> 3] |= 1 << (j & 7);
            if (++mismatches > 0 /*mismatches_allowed*/)
                break;
        }
        ref_nt = p-- & 3 ? ref_nt >> 2 : *--s;
    }
    return mismatches;
}

// returns no mismatches
static int
match_template(gzin_t* gz, struct map_t &map)
{
    // There can still be mismatches in upseq or post key.
    unsigned mismatches = 0;

    uint32_t p = (map.p >> 1) - map.iend;
    uint8_t* s = map.s + (p >> 2);
    char* upseq = map.upseq;
    unsigned ref_nt = *s >> ((p & 3) << 1);

    for (unsigned j = 0; j != map.readlength; ++j) {

        if (j >= map.iend) {
            int c = GZTC(gz);
            if (c == '\n' || c == -1)
                return -1;
            *upseq = toupper(c);
        }

        unsigned read_nt = 0;
        switch(*upseq++) {
            case 'C': read_nt = 2;
            case 'G': read_nt |= 1;
            case 'U':
            case 'T': read_nt ^= 2;
        }

        if (read_nt != (ref_nt & 3)) {

            // TODO: mism[j >> 3] |= 1 << (j & 7);
            if (++mismatches > 0 /*mismatches_allowed*/)
                break;
        }
        ref_nt = ++p & 3 ? ref_nt >> 2 : *++s;
    }
    return mismatches;
}

static int
writeseq(Key_t* kc, gzin_t* gz, struct map_t &map)
{

    char unknown[] = "*";
    char* mtd = unknown;
    char* tid = unknown;
    char* upseq = map.upseq;
    unsigned tln = 0, mps = 0, mq = 0, flag = 44;
    if (map.p) {
        mq = 37;
        flag = map.p & 1;
        flag = (flag << 4) | 0x42;
        tid = kc->id + kc->h[map.ho].ido;
        // header was already printed.

        // XXX: add .corr to map.p ! (after sorting?)
    }

    OPR0("%u\t%s\t%u\t%u\t%uM\t%s\t%u\t%d\t",
        flag, tid, map.p >> 1, mq, map.readlength, mtd, mps, tln);

    // if the orientation is the other, revcmp sequence and reverse qual.
    if (map.p & 1) {

        upseq += map.readlength;// Revcmp orientation
        for (unsigned j = 0; j != map.readlength; ++j) {
            int c = *--upseq;
            if (c == 'A' || c == 'T')
                c ^= 'A' ^ 'T';
            else if (c == 'C' || c == 'G')
                c ^= 'C' ^ 'G';
            fputc(c, stdout);

            c = GZTC(gz);
            if (c == '\n' || c == -1)
                return -1;
            *upseq = c;
        }
        OPR("\t%s\tMD:Z:%u\tNM:i:0", upseq, map.readlength);
    } else {
        OPR0("%s\t", upseq);// Template orientation
        for (unsigned j = 0; j != map.readlength; ++j) {
            int c = GZTC(gz);
            if (c == '\n' || c == -1)
                return -1;
            fputc(c, stdout);
        }
        OPR("\tMD:Z:%u\tNM:i:0", map.readlength);
    }
    return 0;
}


// find minimal key, hopefully unique within extension.
static unsigned
getmink(Key_t* kc, gzin_t* gz, struct map_t &map)
{
    uint32_t* hkoffs = kc->hkoffs + kc->hkoffs_l - 1;
    char* upseq = map.upseq;
    Hdr* h = kc->h + kc->h_l - 1;
    unsigned i = 0, iend = map.readlength;
    int extension = map.readlength - KEY_WIDTH;
    uint32_t t, dna = 0u, rc = 0u;

    uint32_t prevko = -1u;
    const uint64_t last_hs = (kc->s_l >> 2) - h->len;
    map.s = kc->s + last_hs;

    print_hdr(gz);
    EPR("Contig: %s, extension: %u (init)", kc->id + h->ido, extension);
    while (i != iend) {
        int c = GZTC(gz);
        if (c == '\n' || c == -1)
            return -1u;

        // need to preserve N's (if really necessary could use quals though)
        *upseq++ = c = toupper(c);

        NB(c != -1);
        rc = (rc << 2) & KEYNT_MASK;
        dna >>= 2;
        t = 0;
        switch(c) {
            case 'C': t = 2;
            case 'G': t |= 1;
            case 'U':
            case 'T': t ^= 2;
            case 'A':
            default: rc ^= t;
                dna ^= t << KEYNT_TOP;

                if (++i < KEY_WIDTH) // first complete key
                    continue;

                rc ^= dna;            /* rc becomes all deviant bits */
                t = rc & -rc;         /* isolate first deviant bit (devbit) */
                t |= !t;              /* for palindromes use first bit. have to use one.*/
                t = !(t & dna);       /* was devbit not set in dna? */
                // t is orientation of key, 1 if dna, not rc, has same orientation as key

                uint32_t nx = dna ^ (rc & -t); /* dna or rc dependent on devbit */
                rc ^= dna;            /* restore rc */
                nx ^= (-!!(nx & KEYNT_BUFSZ)) & SNDX_TRUNC_MASK; /*shorten index by one */

                NB(nx < KEYNT_BUFSZ);
                uint32_t ko = kc->contxt_idx[nx];

                if (ko >= kc->kct_l) {
                    // non-existent key. There must be a mimsmatch in this key
                    // TODO: use extent above kc->kct_l to indicate possible mismatches
//                    unkey[i >> 3] |= 1 << (i & 7); (or set bit in upseq?)
                    continue;
                }
                if (ko >= prevko)
                    continue; // not a k-mer minimum

                if (ko < kc->hkoffs[kc->h_l])
                    continue; // multimapper k-mer;

                prevko = ko;

                // if t is set, key orientation corresponds with read orientation
                // if kc->kct[ko] & 1 is set, key orientation is template orientation
                // so if map.p & 1 == 0, read orientation is template orientation
                map.p = ((t & 1) ^ (kc->kct[ko] & ~DUP_BIT));

                // if the DUP_BIT is set this indicates a unique key was excised; there
                // should then be another key on this read, also indicating this position.

                // hkoffs indicates how many k's per extension per contig.
                // TODO: rather than per contig, iterate per u32 for multiple contigs.
                EPR("hit");
                while (ko + 1 < *hkoffs) {
                    EPR("ko: %u, *hkoffs:%u", ko, *hkoffs);
                    NB(hkoffs - kc->hkoffs);
                    --hkoffs;
                    if (h - kc->h) {
                        map.s -= (--h)->len;
                    } else {
                        --extension;
                        h = kc->h + kc->h_l - 1;
                        map.s = kc->s + last_hs;
                    }
                    EPR("Contig: %s, extension: %u", kc->id + h->ido, extension);
                }
                EPR("end:ko: %u, *hkoffs:%u", ko, *hkoffs);
                // potential hit;
                EPR("%u, i:%u, extension:%u", iend, i, extension);
                iend = i + extension;
        }
    }
    if (iend < map.readlength) {
        unsigned mismatches;
        map.ho = (hkoffs - kc->hkoffs) % kc->h_l;
        map.iend = iend;
        map.p -= NT_WIDTH;
        EPR(">%s:%u", kc->id + kc->h[map.ho].ido, map.p >> 1);
        // possible hit

        // for orientation, see notes in getmink()
        if (map.p & 1)
            mismatches = match_revcmp(gz, map);
        else
            mismatches = match_template(gz, map);

        skip_plus_line(gz);

        if (mismatches) {
            EPR("got mismatches");

            NB(mismatches != -1u, "sequence truncated?");

            // TODO: try to resolve mismatches.. only if it fails:
            map.p = 0;
        }

    } else { // XXX XXX we end up here?
        EPR("reached end.");
        NB (iend != -1u, "truncated sequence?");
        map.p = 0;
    }
    return writeseq(kc, gz, map);
}

static int
fq_read(Key_t* kc, gzin_t* gz)
{
    map_t map = {0};
    int res, c = map.readlength = kc->readlength;

    // upseq: to store sequence before key + null.
    map.upseq = (char*)malloc(c + 1);
    map.upseq[c] = '\0';

    c = (c >> 3) + !!(c & 3);
    /*map.mism = (uint8_t*)malloc(c);
    map.unkey = (uint8_t*)malloc(c);

    while (c--)
        map.unkey[c] = map.mism[c] = 0;*/

    //struct mapstat_t ms = {0};
    do {
        _EVAL(getmink(kc, gz, map));

    } while (GZTC(gz) == '\n');
    c = (c == -1); // end of file
out:
    QARN(c < 0, "error:%d", c);
    //free(map.mism);
    //free(map.unkey);
    free(map.upseq);
    return c;
}

int
map_fq_se(struct seqb2_t* sb2, char C*C cmdl)
{
    int res = -ENOMEM;

    struct gzfh_t* fhio[2] = { sb2->fh + 1, sb2->fh + 2};
    const char* ext[2] = {".uq",".2b"};
    char file[512];
    Key_t kc = {0};
    gzin_t gz = {0};

    kc.readlength = sb2->readlength;
    ASSERT(fhio[1]->name != NULL, return -EFAULT);
    unsigned len = strlen(fhio[1]->name) + 1;
    ASSERT(strstr(fhio[1]->name, ext[1]), return -EFAULT);
    // TODO: first read in several reads to verify

    _ACTION(init_fq(sb2), "intializing memory");

    for (int i=0; i != 2; ++i) {
        if (fhio[i]->name == NULL) {
            fhio[i]->name = &file[len*i];
            strncpy(fhio[i]->name, fhio[1]->name, len);
            _ACTION0(reopen(fhio[i], ext[1], ext[i]),
                    "%s %s", ext[i], fhio[i]->fp ? "exists" : "does not exist");
        }
    }
    ASSERT(fhio[0]->fp && fhio[1]->fp, return -EFAULT,
            "need seqb2, keycount and unique boundary files");
    kc.contxt_idx = buf_init_arr(kc.contxt_idx, KEYNT_BUFSZ_SHFT);
    // 2) open seqb2 for verification of reads
    // 3) open kc
    _ACTION(load_seqb2(fhio[1], &kc), "loading twobit sequence file");
    _ACTION(load_kc(fhio[0], &kc), "loading keycounts file");

    // 4) print header
    print_hdr_lines(&kc, cmdl);
    // 5) open fq for reading
    set_readfunc(sb2->fh, &gz);
    _ACTION(fq_read(&kc, &gz), "mapping reads");
    //...
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free(sb2->lookup);
    free(sb2->s);
    free_kc(&kc);
    return res;
}
