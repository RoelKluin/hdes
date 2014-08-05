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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <getopt.h>
#include "seq.h"
#include "fq.h"
#include "gz.h"
#include "util.h"

static int
b2_write(gzfh_t *fh, char *start, int len)
{
    if ((len = write(fileno(fh->fp), start, len)) < 0) {
        fputs("error while writing\n", stderr);
        return -3;
    }
    if (DEBUG)
        fprintf(stderr, "%u bytes written\n", len);

    return ferror(fh->fp) ? -3 : len;
}

typedef struct option_description
{
    const char *name;
    int has_arg;
    int *flag;
    int val;
    const char *descr;
} option_description_t;

// order matters and should first be input, then output, then integer args, then etc.
static const char optstr[] = "1:2:r:o:m:l:b:p:fh";
static struct option_description dopt[] = {
    {"fastq", required_argument, NULL, '1', "<FILE>\t[first] input fastq"},
    {"2nd-fastq", required_argument, NULL, '2', "<FILE>\t2nd input fastq (optional)"},
    {"reference", required_argument, NULL, 'r', "<FILE>\tinput fasta reference (optional)"},
    {"out", required_argument, NULL, 'o', "<FILE>\toutput file prefix"},
    {"readlength", required_argument, NULL, 'l', "<INT>\treadlength (optional)"},
    {"maxreads", required_argument, NULL, 'm', "<INT>\tonly this many reads (optional)"},
    {"blocksize", required_argument, NULL, 'b', "<INT>\tgzip blocksize in Mb"},
    {"phred", required_argument, NULL, 'p', "<INT>\tphred offset for quality scores"},
    {"force", no_argument, NULL, 'f', "\tforce overwrite"},
    {"help", no_argument, NULL, 'h', "\tprint options"},
    {NULL, 0, NULL, 0, NULL}
};

static int usage()
{
    unsigned i;
    fprintf(stderr, "Program: %s\n", PROGRAM_NAME);
    fprintf(stderr, "Version: %s\n", PROGRAM_VERSION);
    fprintf(stderr, "Contact: Roel Kluin <r.kluin@nki.nl>\n\n");
    fprintf(stderr, "Usage:   %s [options] <in.fq> [in2.fq] [ref] <prefix> \n\n",
            PROGRAM_NAME);
    for (i = 0; dopt[i].name != NULL; ++i)
        fprintf(stderr, "\t-%c|--%s\t%s\n", dopt[i].val, dopt[i].name, dopt[i].descr);
    fputc('\n', stderr);
    return 1;
}

int main(int argc, char* const* argv)
{
    struct seqb2_t seq = {0}; /* init everything to 0 or NULL */
    unsigned t, i = 0u, fhsz = ARRAY_SIZE(seq.fh);
    int c, ret = EXIT_FAILURE;
    int optvals[][4] = {
        // (default) value, bottom limit, upper lim, enforce
        {'l', 0, KEY_WIDTH, SEQ_OFFSET_MAX},
        {'m', -1, 1, 0x7fffffff},
        {'b', 32, 1, 8192},
        {'p', 33, 32, 89}
    };
    seq.fh[0].fp = seq.fh[2].fp = stdin;
    seq.fh[3].fp = stdout;
    seq.fh[3].write = &b2_write;

    fputc('\n', stderr);

    /* parse cmdline args */
    while ((c = getopt_long_only(argc, argv, optstr, (option*)dopt, (int*)&i)) != -1) {
        i ^= i;
        switch (c) {
            case 'o': ++i;
            case 'r': ++i;
            case '2': ++i;
            case '1':
                fprintf(stderr, "%s: %s\n", dopt[i].name, optarg);
                assert((seq.fh[i].write == NULL) != (i == fhsz - 1));
                seq.fh[i].name = optarg;
                break;
            case 'p': ++i;
            case 'b': ++i;
            case 'm': ++i;
            case 'l': assert(dopt[i + fhsz].val == optvals[i][0]);
//                fprintf(stderr, "%s: %s\n", dopt[i + fhsz].name, optarg);
                c = atoi(optarg);
                t = (c < optvals[i][2]);
                if (t || (c > optvals[i][3])) {
                    t = i + fhsz;
                    fprintf(stderr, "\t-%c|--%s %s is invalid.\n\t%s must be "
                        "between %u and %u\n", dopt[t].val, dopt[t].name, optarg,
                        dopt[t].descr, optvals[i][2], optvals[i][3]);
                    goto out;
                }
                optvals[i][1] = c;
                break;
            case 'h': usage();
            case '?': case ':': fputc('\n', stderr); goto out;
            default:  seq.mode |= amopt(c);
        }
    }
    seq.readlength = optvals[0][1];
    seq.readlimit = optvals[1][1];
    seq.phred_offset = optvals[3][1];

    /* open files and allocate memory */
    for (i=0; i != fhsz; ++i)
        if (seq.fh[i].name == NULL && optind != argc)
            seq.fh[i].name = argv[optind++];

    if (optind != argc) {
        while (optind != argc)
            fprintf(stderr, "main: Extra argument %s\n", argv[optind++]);
        fprintf(stderr, "main: Too many arguments\n");
        goto out;
    }

    for (i=0; i != fhsz; ++i) {
        if (seq.fh[i].fp == NULL) {
            if ((dopt[i].val == '2') || (dopt[i].val == 'r')) continue; // SE or assembly

            fprintf(stderr, "main: Missing -%c [%s].\n",
                    dopt[i].val, dopt[i].name);
            goto out;
        }
        if (set_io_fh(&seq.fh[i], seq.mode, (uint64_t)optvals[0][0] << 20) != 0) {
            fprintf(stderr, "main: -%c [%s] failed.\n",
                    dopt[i].val, dopt[i].name);
            goto out;
        }
    }

    if ((c = init_fq(&seq)) < 0) {
        fprintf(stderr, "ERROR: init_fq() returned %d\n", c);
        goto out;
    }
    if ((c = fq_b2(&seq)) != 0) {
        fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
        goto out;
    }
    fq_print(&seq);
    /* read fasta, transform and write file */
/*    c = fa_b2(&seq, readlength);
    if (c < 0)
        fprintf(stderr, "main: fa_b2() returned %d\n", c);
*/
    /*c = seq_b2(&seq);
    if (c < 0)
        fprintf(stderr, "main: fa_b2() returned %d\n", c);*/
    ret = EXIT_SUCCESS;
out: /* cleanup */
    free_fq(&seq);
    for (i=0; i != fhsz; ++i) {
        if (seq.fh[i].io == NULL) continue;
        // XXX: valgrind complains here but the problem is in zlib
        // probably not a bug.
        if (seq.fh[i].close && seq.fh[i].close(seq.fh[i].io) != Z_OK) {
            fprintf(stderr, "main: gzclose fails for %s\n", dopt[i].name);
            close(seq.fh[i].fd);
        }
        close(seq.fh[i].fd);
    }
    return ret;
}

/*int main(int argc, const char* argv[])
{
    using namespace std;
    int ret = EXIT_FAILURE;
    igzstream fa(argv[1]);
    igzstream fq(argv[2]);
    unsigned phred_offset = 33;

    if (argc > 4) {
        phred_offset = atoi(argv[3]);
        cerr << "Phred offset set to " << phred_offset << endl;
    }
    assert((phred_offset >= 33) && (phred_offset <= 126));
    unsigned readlength = atoi(argv[argc - 1]);
    struct fq_arr uqct(readlength, phred_offset);

    string sq;

    while (getline(fa, sq)) {
        char c = sq[0];
        if (c == '>') {
            // refill, dna & revcmp up to KEY_WIDTH Nts
            uqct.pre_init_key();
            while(uqct.init_key(seq) == false)
                if (not getline(fa, sq)) goto err;

            uqct.add_max_key(sq + KEY_WIDTH);
       }
        uqct.fa_put((const uint8_t*)sq.c_str());
        if (uqct.is_err) goto err;
    }

    string nm, tmp, ql;
    while (getline(fq, nm) && getline(fq, sq) && getline(fq, tmp) && getline(fq, ql)) {
        uqct.put((const uint8_t*)sq.c_str(), (const uint8_t*)ql.c_str(), (const uint8_t*)nm.c_str());
        if (uqct.is_err) goto err;
    }

err: uqct.dump();
    fq.close();
    fa.close();
    return ret;
}*/


