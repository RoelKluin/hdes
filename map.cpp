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

/*#define GZTC(gz) ({\
  int __c = ((gz)->c((gz)->g));\
  if (__c != -1)\
      EPR0("%c", (char)__c);\
  else\
    EPR("<EOF>");\
  __c;\
})*/

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
        OPR("@SQ\tSN:%s\tLN:%u", kc->id + h->ido, h->end);
        EPR("@SQ\tSN:%s\tLN:%u", kc->id + h->ido, h->end);
    }
    OPR("@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s", PROGRAM_NAME, PROGRAM_NAME,
            PROGRAM_VERSION, commandline);
}

static inline void
print_hdr(gzin_t* gz)
{
    int c;
    while ((c = GZTC(gz)) != '\n') {
        NB(c != -1);
        fputc(c, stdout);
    }
    fputc('\t', stdout);
}

static inline unsigned
match_nt(char rdchar, unsigned ref_nt, unsigned read_nt)
{
    switch(rdchar) {
        case 'C': read_nt ^= 2;
        case 'G': read_nt ^= 1;
        case 'U':
        case 'T': read_nt ^= 2;
    }
    //Alignment;
    //EPR("%c <-> %c", (char)b6(rdb2 << 1), (char)b6((refb2&3) << 1));

    return read_nt != (ref_nt & 3);
}

// Als de read korter is dan map.readlength dan wordt dit hier gecorrigeerd.
static void
match_revcmp(struct map_t &map)
{
    NB(map.p >= (KEY_WIDTH + map.readlength - map.i) << 1);
    // There can still be mismatches in upseq or post key.
    map.p += (map.i - KEY_WIDTH) << 1;
    uint32_t p = map.p >> 1; // set to last 2bit.

    uint8_t* s = map.s + (p >> 2);
    char* upseq = map.upseq;
    unsigned ref_nt = *s >> ((p & 3) << 1);
    unsigned j;

    for(j = 0; j != map.readlength; ++j) {

        //EPR0("p:%u(%u)\tso:%lu\t", p, j, s - map.s);
        unsigned mm = match_nt(*upseq++, ref_nt, 2);

        map.mismatches += mm;
        // TODO: mism[j >> 3] |= mm << (j & 7); (or set bits in upseq)

	ref_nt = --p & 3 ? ref_nt >> 2 : *--s;
    }
    *upseq = '\0';
    map.p -= (map.readlength - 1) << 1; // SAM is one-based.
}

// returns no mismatches
static void
match_template(struct map_t &map)
{
    // There can still be mismatches in upseq or post key.

    NB(map.p >= (map.i << 1));
    map.p -= (map.i << 1);
    uint32_t p = map.p >> 1;
    map.p += 2; // SAM is one-based, and map.p.
    uint8_t* s = map.s + (p >> 2);

    char* upseq = map.upseq;
    unsigned ref_nt = *s >> ((p & 3) << 1);
    unsigned j;

    for(j = 0; j != map.readlength; ++j) {

        unsigned mm = match_nt(*upseq++, ref_nt, 0);

        map.mismatches += mm;
        // TODO: mism[j >> 3] |= mm << (j & 7); (or set bits in upseq)

        ref_nt = ++p & 3 ? ref_nt >> 2 : *++s;
    }
    *upseq = '\0';
}


static void
writeseq(Key_t* kc, gzin_t* gz, struct map_t &map)
{

    char unknown[] = "*";
    char* mtd = unknown;
    char* tid = unknown;
    char* upseq = map.upseq;
    unsigned tln = 0, mps = 0, mq = 0, flag = 44;
    if (map.mismatches != -1) {
        mq = 37;
        flag = map.p & 1;
        flag = (flag << 4) | 0x42;
        tid = kc->id + kc->h[map.ho].ido;
        // header was already printed.

        // XXX: add .corr to map.p ! (after sorting?)
    }

    // skip plus line
    int c;
    NB(GZTC(gz) == '+');
    while ((c = GZTC(gz)) != '\n')
        NB(c != -1);

    OPR0("%u\t%s\t%u\t%u\t%uM\t%s\t%u\t%d\t",
        flag, tid, map.p >> 1, mq, map.readlength, mtd, mps, tln);

    // if the orientation is the other, revcmp sequence and reverse qual.
    if (map.p & 1) {
        //EPR("// Revcmp orientation");
        upseq += map.readlength;
        for (unsigned j = 0; j != map.readlength; ++j) {
            c = *--upseq; // & 0x5f discard mismatch and non-existent key markings.. template?
            if (c == 'A' || c == 'T')
                c ^= 'A' ^ 'T';
            else if (c == 'C' || c == 'G')
                c ^= 'C' ^ 'G';
            fputc(c, stdout);

            c = GZTC(gz);
            NB(c != '\n', "Quality line too short (rc).");
            NB(c != -1, "Fastq truncated in quality (rc).");
            *upseq = c; // replace by qualities
        }
        OPR("\t%s\tMD:Z:%u\tNM:i:%u\tXE:i:%u\tXI:i:%u\tXN:i:%u", upseq, map.readlength, map.mismatches, map.extension, map.iteration, map.nonk);
    } else {
        //EPR("// Template orientation");
        OPR0("%s\t", upseq);
        for (unsigned j = 0; j != map.readlength; ++j) {
            c = GZTC(gz);
            NB(c != '\n', "Quality line too short.");
            NB(c != -1, "Fastq truncated in quality.");
            fputc(c, stdout);
        }
        OPR("\tMD:Z:%u\tNM:i:%u\tXE:i:%u\tXI:i:%u\tXN:i:%u", map.readlength, map.mismatches, map.extension, map.iteration, map.nonk);
    }
}


// find minimal key, hopefully unique within extension.
static unsigned
getmink(Key_t* kc, gzin_t* gz, struct map_t &map)
{
    uint32_t* hkoffs = kc->hkoffs + kc->hkoffs_l - 1;
    char* upseq = map.upseq;
    Hdr* h = kc->h + kc->h_l - 1;
    unsigned i = 0, iend = -1u;

    uint32_t* ext_iter = kc->ext_iter + kc->ext_iter_l - 1;
    unsigned iteration = *ext_iter;

    uint32_t t, dna = 0u, rc = 0u;

    uint32_t prevko = -1u;
    const uint64_t last_hs = (kc->s_l >> 2) - h->len;
    map.s = kc->s + last_hs;

    int c = GZTC(gz);
    if (c == -1) return 0; // end of fastq
    NB(c == '@', "Expected fastq header '@', not '%c'", (char)c);
    print_hdr(gz);

//    EPR("Contig: %s, extension: %u, iteration: %u (init)", kc->id + h->ido, extension, iteration);
    while ((c = GZTC(gz)) != '\n') {

        NB (c != -1, "Fastq truncated in sequence?");

        // need to preserve N's (if really necessary could use quals though)
        *upseq++ = c = toupper(c);
        if (++i >= iend) continue; // altijd seq laten aflopen tot aan '\n' XXX >= ?

        rc = (rc << 2) & KEYNT_MASK;
        dna >>= 2;
        t = 0;
        switch(c) {
            case 'C': t ^= 2;
            case 'G': t ^= 1;
            case 'U':
            case 'T': t ^= 2;
            default: rc ^= t ^ 2;
                dna ^= t << KEYNT_TOP;
	}

	if (i < KEY_WIDTH) // first complete key
	    continue;

	rc ^= dna;            /* rc becomes all deviant bits */
	t = rc & -rc;         /* isolate first deviant bit (devbit) */
	t |= !t;              /* for palindromes use first bit. have to use one.*/
	t = !(t & dna);       /* was devbit not set in dna? */
	// t is orientation of key, 1 if dna, not rc, has same orientation as key

	uint32_t nx = dna ^ (rc & -t); /* dna or rc dependent on devbit */
	rc ^= dna;            /* restore rc */
	nx ^= (-!!(nx & KEYNT_BUFSZ)) & SNDX_TRUNC_MASK; /*shorten index by one */

	NB(nx < KEYNT_BUFSZ, "nx too long");
	uint32_t ko = kc->contxt_idx[nx];

	if (ko >= kc->kct_l) {
	    // non-existent key. There must be a mimsmatch in this key
	    // TODO: use extent above kc->kct_l to indicate possible mismatches
//                    unkey[i >> 3] |= 1 << (i & 7); (or set bit in upseq?)
	    map.nonk++;
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
	map.p = (t & 1) ^ (kc->kct[ko] & ~DUP_BIT);

	// if the DUP_BIT is set this indicates a unique key was excised; there
	// should then be another key on this read, also indicating this position.

	// hkoffs indicates how many k's per extension per contig.
	// TODO: rather than per contig, iterate per u32 for multiple contigs.
	if (ko < *hkoffs) {
	    while (ko < *(hkoffs - kc->h_l - 1)) {

		if (iteration-- == 0)
		    iteration = *--ext_iter;

		//R;
		hkoffs -= kc->h_l;
		if (hkoffs - kc->hkoffs < kc->h_l)
		    break;
	    }
	    while (ko < *(hkoffs - 1)) {
		NB(hkoffs - kc->hkoffs, "hkoffs out of range");
		--hkoffs;
		//EPR("ko: %u, *hkoffs:%u", ko, *hkoffs);
		if (h - kc->h) {
		    --h;//R;
		    map.s -= h->len;
		} else {
		    if (iteration-- == 0)
			iteration = *--ext_iter;

		    h = kc->h + kc->h_l - 1;
		    map.s = kc->s + last_hs;
		}
//                    EPR("Contig: %s, extension: %u, iteration: %u", kc->id + h->ido, extension, iteration);
	    }
	}
	//for (uint32_t* z = hkoffs; z - kc->hkoffs != -1; --z)
	//    EPR("hkoffs[%lu] = %u", z - kc->hkoffs, *z);
	//R;
//                EPR("end:ko: %u, *hkoffs:%u", ko, *hkoffs);
	// potential hit;
//                EPR("%u, i:%u, extension:%u, iteration: %u", iend, i, extension, iteration);
	map.extension = ext_iter - kc->ext_iter;
	map.iteration = iteration;
	map.i = i;
	iend = i + (ext_iter - kc->ext_iter - 1);
    }
    map.readlength = i;
    map.iend = iend;
    map.mismatches = 0;

    if (iend <= map.readlength) {
        map.ho = (hkoffs - kc->hkoffs) % kc->h_l;
        map.iend = iend;
//        EPR(">%s:%u", kc->id + kc->h[map.ho].ido, map.p >> 1);
        // possible hit

        // for orientation, see notes in getmink()
        if (map.p & 1)
            match_revcmp(map);
        else
            match_template(map);

        //if (map.mismatches > 0 /*mismatches_allowed*/) {

            // TODO: try to resolve mismatches.. only if it fails:
        //    map.p = 0; // Mismatches
        //}

    } else {
	map.mismatches = -1;//multimapper

        NB (map.readlength != 0, "Expected sequence, but the line was empty.");
    }
    writeseq(kc, gz, map);
    return 1;
}

static int
fq_read(Key_t* kc, gzin_t* gz)
{
    map_t map = {0};
    int c = map.readlength = kc->readlength;

    // upseq: to store sequence before key + null.
    map.upseq = (char*)malloc(c + 1);
    map.upseq[c] = '\0';

    /*c = (c >> 3) + !!(c & 3);
    map.mism = (uint8_t*)malloc(c);
    map.unkey = (uint8_t*)malloc(c);

    while (c--)
        map.unkey[c] = map.mism[c] = 0;*/

    //struct mapstat_t ms = {0};
    while (getmink(kc, gz, map)) {
        c = GZTC(gz);
        if (c == -1)
            break;
        NB(c == '\n');
    }
    c = 0; // end of file

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
