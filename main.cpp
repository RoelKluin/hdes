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
#include <getopt.h>
#include <string.h> //strstr()
#include "fa.h"
#include "fq.h"

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
    {"ref", required_argument, NULL, 'r', "<FILE>\tinput fasta reference (optional)"},
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
    fputs("Instead of a file you may use stdin or stdout\n", stderr);
    fputc('\n', stderr);
    return 1;
}

int main(int argc, char* const* argv)
{
    struct seqb2_t seq = {0}; /* init everything to 0 or NULL */
    unsigned t, i = 0u, fhsz = ARRAY_SIZE(seq.fh);
    int c, ret = EXIT_FAILURE;
    int optvals[][4] = {
        // (default) value, bottom limit, upper lim, option
        {0, KEY_WIDTH, SEQ_OFFSET_MAX, 'l'},
        {-1, 1, 0x7fffffff, 'm'},
        {32, 1, 1024, 'b'},
        {33, 32, 89, 'p'}
    };
//    seq.fh[0].fp = seq.fh[2].fp = stdin; // default to stdin
    seq.fh[fhsz - 1].fp = stdout;
    seq.fh[fhsz - 1].write = &b2_write;


    /* parse cmdline args */
    while ((c = getopt_long_only(argc, argv, optstr, (option*)dopt, (int*)&i)) != -1) {
        i ^= i;
        switch (c) {
            case 'o': ++i;
            case 'r': ++i;
            case '2': ++i;
            case '1':
//                fprintf(stderr, "%s: %s\n", dopt[i].name, optarg);
                assert((seq.fh[i].write == NULL) != (i == fhsz - 1));
                seq.fh[i].name = optarg;
                break;
            case 'p': ++i;
            case 'b': ++i;
            case 'm': ++i;
            case 'l': assert(dopt[i + fhsz].val == optvals[i][3]);
//                fprintf(stderr, "%s: %s\n", dopt[i + fhsz].name, optarg);
                c = atoi(optarg);
                t = (c < optvals[i][1]);
                if (t || (c > optvals[i][2])) {
                    t = i + fhsz;
                    fprintf(stderr, "\t-%c|--%s %s is invalid.\n\t%s must be "
                        "between %u and %u\n", dopt[t].val, dopt[t].name, optarg,
                        dopt[t].descr, optvals[i][1], optvals[i][2]);
                    goto out;
                }
                optvals[i][0] = c;
                break;
            case 'h': usage();
            case '?': case ':': fputc('\n', stderr); goto out;
            default:  seq.mode |= amopt(c);
        }
    }
    seq.readlength = optvals[0][0];
    seq.readlimit = optvals[1][0];
    seq.blocksize = optvals[2][0];
    seq.phred_offset = optvals[3][0];

    /* options without a flag */
    while (optind != argc) {
        char* f = (char*)argv[optind++];
        i = get_fastx_type(f, 2, fhsz);
        while ((i < 2) && seq.fh[i].name != NULL) ++i; // get available read fastqs, skipped for fasta

        if ((i >= fhsz) || seq.fh[i].name != NULL) { // fhsz , not neccearily one bit, marks error.
            fprintf(stderr, "Unhandled argument %s (%u)\n", f, i);
            goto out;
        }
        seq.fh[i].name = f;
    }

    if (optind != argc) {
        while (optind != argc)
            fprintf(stderr, "main: Extra argument %s\n", argv[optind++]);
        fprintf(stderr, "main: Too many arguments\n");
        goto out;
    }

    /* open files and allocate memory */
    for (i=0; i != fhsz; ++i) {
        if (seq.fh[i].name == NULL) continue;
        c = set_io_fh(&seq.fh[i], &seq.mode, (uint64_t)seq.blocksize << 20);
        if (c < 0) {
            fprintf(stderr, "main: -%c [%s] failed.\n", dopt[i].val, dopt[i].name);
            goto out;
        }

        fprintf(stderr, "%s:\t", dopt[i].name);
        if (c == 0)
            fprintf(stderr, "%s\n", seq.fh[i].name);
        else if (c == 1) fputs("(none given)\n", stderr);

    }

    if (seq.fh[2].name != NULL) {
        if (seq.fh[1].name) {
            fputs("Paired-end alignment\n", stderr);
            /*if ((c = pe_fq_b2(&seq)) != 0) {
                fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);*/
        } else if (seq.fh[0].name) {
            fputs("Single-end alignment\n", stderr);
            /*if ((c = fq_b2(&seq)) != 0) {
                fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);*/
        } else {
            //if (seq.readlength == 0) {
            //    fputs("Readlength is required for indexing\n", stderr);
            //    goto out;
            //}
            fputs("Indexing\n", stderr);
            if ((c = fa_ndx(&seq)) != 0) {
                fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
                goto out;
            }
            fputc('\n', stderr);
        }
    } else {
        if (seq.fh[1].name) {
            fputs("Paired-end assembly\n", stderr);
            /*if ((c = pe_fq_b2(&seq)) != 0) {
                fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);*/
        } else if (seq.fh[0].name) {
            fputs("Single-end assembly\n", stderr);
            if ((c = fq_b2(&seq)) != 0) {
                fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);
        } else {
            fputs("No input files specified\n", stderr);
            goto out;
        }
    }

    ret = EXIT_SUCCESS;
out: /* cleanup */
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


