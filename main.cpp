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

#ifndef PROGRAM_NAME
# define PROGRAM_NAME "uqct"
#endif
#ifndef PROGRAM_VERSION
# define PROGRAM_VERSION "0.009"
#endif

static int
b2_write(const gzfh_t *fh, const char *s, uint64_t l,
        const size_t sz)
{
    while (l) {
        int c = write(fileno(fh->fp), s, l > INT_MAX ? INT_MAX : l);
        if (c < 0) {
            fputs("error while writing\n", stderr);
            return c;
        }
        fprintf(stderr, "==%d bytes written\n", c);
        c /= sz, l -= c, s += c;
    }
    return ferror(fh->fp) ? -3 : 0;
}

typedef struct option_description
{
    const char *name;
    int has_arg;
    int *flag;
    unsigned val;
    const char *descr;
} option_description_t;

// order matters and should first be input, then output, then integer args, then etc.
static const char optstr[] = "1:2:r:o:m:l:b:p:fh";
static struct option_description dopt[] = {
    {"fastq", required_argument, NULL, '1', "<FILE>\t[first] input fastq"},
    {"2nd-fastq", required_argument, NULL, '2', "<FILE>\t2nd input fastq (optional)"},
    {"ref", required_argument, NULL, 'r', "<FILE>\tinput fasta reference (optional)"},
    {"out", required_argument, NULL, 'o', "<FILE>\toutput file prefix"},
    {"maxreads", required_argument, NULL, 'm', "<INT>\tonly this many reads (optional)"},
    {"readlength", required_argument, NULL, 'l', "<INT>\treadlength (optional)"},
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
    fputs("Instead of a file you may specify stdin or stdout for some tasks\n", stderr);
    fputc('\n', stderr);
    return 1;
}

int main(int argc, char* const* argv)
{
    uint64_t blocksize;
    struct seqb2_t seq = {0}; /* init everything to 0 or NULL */
    unsigned i = 0u, fhsz = ARRAY_SIZE(seq.fh);
    int c, ret = EXIT_FAILURE;
    uint32_t optvals[16] = {
        // maxreads, readlength, blocksize, phred_offset
        -1u, 0, 32, 33,                // default values
        1, KEY_WIDTH, 1, 32,           // bottom limit (minima)
        -1u, SEQ_OFFSET_MAX, 1024, 89, // maxima
        'm', 'l', 'b', 'p'             // respective option
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
            case 'l': ++i;
            case 'm': {
                assert(dopt[i + fhsz].val == optvals[12 + i]);
//                fprintf(stderr, "%s: %s\n", dopt[i + fhsz].name, optarg);
                unsigned t = atoi(optarg);
                if ((t < optvals[4 + i]) || (t > optvals[8 + i])) {
                    t = i + fhsz;
                    fprintf(stderr, "\t-%c|--%s %s is invalid.\n\t%s must be "
                        "between %u and %u\n", dopt[t].val, dopt[t].name, optarg,
                        dopt[t].descr, optvals[4 + i], optvals[8 + i]);
                    goto out;
                }
                optvals[i] = t;
                break;
            }
            case 'h': usage();
            case '?': case ':': fputc('\n', stderr); goto out;
            default:  seq.mode |= amopt(c);
        }
    }
    seq.maxreads = optvals[0];
    seq.readlength = optvals[1];
    seq.blocksize = optvals[2];
    seq.phred_offset = optvals[3];
    blocksize = (uint64_t)seq.blocksize << 20;

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

        c = set_stdio_fh(&seq.fh[i], &seq.mode);
        if (c == 0)
            c = set_io_fh(&seq.fh[i], blocksize, (seq.mode & amopt('f')) != 0);

        if (c < 0) {
            fprintf(stderr, "main: -%c [%s] failed.\n", dopt[i].val, dopt[i].name);
            goto out;
        }

        fprintf(stderr, "%s(%u):\t%s\n", dopt[i].name, i, c == 0 ? seq.fh[i].name :
                (i >= fhsz - 1 ? "stdout" : "stdin"));
    }
    if (seq.fh[2].name != NULL) {
        if (seq.fh[1].name) {
            fputs("== Paired-end alignment\n", stderr);
            /*if ((c = pe_fq_b2(&seq)) != 0) {
                fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);*/
        } else if (seq.fh[0].name) {
            fputs("== Single-end alignment\n", stderr);
            /*if ((c = fq_b2(&seq)) != 0) {
                fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);*/
        } else {
//            if (seq.readlength == 0) {
//                fputs("== Readlength needed for indexing.\n", stderr);
//                goto out;
//            }
            c = fa_index(&seq);
            if (c < 0) {
                fputs("== failed to create keyct.\n", stderr);
                goto out;
            }
            fputc('\n', stderr);
        }
    } else {
        if (seq.fh[1].name) {
            fputs("== Paired-end assembly\n", stderr);
            /*if ((c = pe_fq_b2(&seq)) != 0) {
                fprintf(stderr, "ERROR: fq_b2() returned %d\n", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);*/
        } else if (seq.fh[0].name) {
            fputs("== Single-end assembly\n", stderr);
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
    fputs("Cleanup\n", stderr); fflush(NULL);

    ret = EXIT_SUCCESS;
out: /* cleanup */
    for (i=0; i != fhsz; ++i) {
        if (seq.fh[i].fp == NULL) continue;
        fprintf(stderr, "closing %u\n", i);
        // XXX: valgrind complains here but the problem is in zlib
        // probably not a bug.
        if (seq.fh[i].close && seq.fh[i].io != Z_NULL && seq.fh[i].close(seq.fh[i].io) != Z_OK)
            fprintf(stderr, "main: gzclose fails for %s\n", dopt[i].name);
        close(seq.fh[i].fd);
    }
    return ret;
}


