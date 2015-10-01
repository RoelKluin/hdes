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
#include <signal.h>
#include "fa.h"
#include "fq.h"

#ifndef PROGRAM_NAME
# define PROGRAM_NAME "uqct"
#endif
#ifndef PROGRAM_VERSION
# define PROGRAM_VERSION "0.016"
#endif

typedef struct option_description
{
    const char *name;
    int has_arg;
    int *flag;
    unsigned val;
    const char *descr;
} option_description_t;

// order matters and should first be input, then output, then integer args, then etc.
static const char optstr[] = "1:2:r:L:o:m:l:b:p:fh";
static struct option_description dopt[] = {
    {"fastq", required_argument, NULL, '1', "<FILE>\t[first] input fastq"},
    {"2nd-fastq", required_argument, NULL, '2', "<FILE>\t2nd input fastq (optional)"},
    {"ref", required_argument, NULL, 'r', "<FILE>\tinput fasta reference (optional)"},
    {"bed", required_argument, NULL, 'L', "<FILE>\ttarget regions bed file (optional)"},
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
    EPR("Program: %s\n"
        "Version: %s\n"
        "Contact: Roel Kluin <r.kluin@nki.nl>\n\n"
        "Usage:   %s [options] [in.fq] [in2.fq] [ref] [bed] [output]\n",
            PROGRAM_NAME, PROGRAM_VERSION, PROGRAM_NAME);
    for (i = 0; dopt[i].name != NULL; ++i)
        EPR("\t-%c|--%s\t%s", dopt[i].val, dopt[i].name, dopt[i].descr);
    EPR("\nInstead of a file you may specify stdin or stdout for some tasks\n");
    return 1;
}

int main(int argc, char* const* argv)
{
    struct seqb2_t seq = {0}; /* init everything to 0 or NULL */
    unsigned i = 0u, fhsz = ARRAY_SIZE(seq.fh);
    int c, ret = EXIT_FAILURE;
    uint32_t blocksize, optvals[16] = {
        // maxreads, readlength, blocksize, phred_offset
        -1u, 0, 32, 33,                // default values
        1, KEY_WIDTH, 1, 32,           // bottom limit (minima)
        -1u, SEQ_OFFSET_MAX, 1024, 89, // maxima
        'm', 'l', 'b', 'p'             // respective option
    };
//    seq.fh[0].fp = seq.fh[2].fp = stdin; // default to stdin
    seq.fh[fhsz - 1].fp = stdout;
    seq.fh[fhsz - 1].write = &b2_write;
    signal(SIGSEGV, handler);   // install our handler


    /* parse cmdline args */
    while ((c = getopt_long_only(argc, argv, optstr, (option*)dopt, (int*)&i)) != -1) {
        i ^= i;
        switch (c) {
            case 'o': ++i;
            case 'L': ++i;
            case 'r': ++i;
            case '2': ++i;
            case '1':
//                EPR("%s: %s", dopt[i].name, optarg);
                assert((seq.fh[i].write == NULL) != (i == fhsz - 1));
                seq.fh[i].name = optarg;
                break;
            case 'p': ++i;
            case 'b': ++i;
            case 'l': ++i;
            case 'm': {
                assert(dopt[i + fhsz].val == optvals[12 + i]);
//                EPR("%s: %s", dopt[i + fhsz].name, optarg);
                unsigned t = atoi(optarg);
                if ((t < optvals[4 + i]) || (t > optvals[8 + i])) {
                    t = i + fhsz;
                    EPR("\t-%c|--%s %s is invalid.\n\t%s must be between %u and %u",
                            dopt[t].val, dopt[t].name, optarg, dopt[t].descr,
                            optvals[4 + i], optvals[8 + i]);
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
    blocksize = optvals[2];
    seq.phred_offset = optvals[3];

    /* options without a flag */
    while (optind != argc) {
        char* f = (char*)argv[optind++];
        i = get_fastx_type(f, fhsz);
        if ((i < 2) && seq.fh[i].name != NULL) ++i; // get available read fastqs, skipped for fasta

        if ((i >= fhsz) || seq.fh[i].name != NULL) { // ge fhsz marks error.
            EPR("Unhandled argument %s (%u)", f, i);
            goto out;
        }
        seq.fh[i].name = f;
    }

    if (optind != argc) {
        while (optind != argc)
            EPR("main: Extra argument %s", argv[optind++]);
        EPR("main: Too many arguments");
        goto out;
    }

    /* open files and allocate memory */
    for (i=0; i != fhsz; ++i) {
        seq.fh[i].blocksize = blocksize;
        if (seq.fh[i].name == NULL) continue;

        c = set_stdio_fh(&seq.fh[i], &seq.mode);
        if (c == 0)
            c = set_io_fh(&seq.fh[i], (seq.mode & amopt('f')) != 0);

        if (c < 0) {
            EPR("main: -%c [%s] failed.", dopt[i].val, dopt[i].name);
            goto out;
        }

        EPR("%s(%u):\t%s", dopt[i].name, i, c == 0 ? seq.fh[i].name :
                (i >= fhsz - 1 ? "stdout" : "stdin"));
    }
    if (seq.fh[2].name != NULL) {
        if (seq.fh[1].name) {
            EPR("== Paired-end alignment");
            /*if ((c = pe_fq_b2(&seq)) != 0) {
                EPR("ERROR: fq_b2() returned %d", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);*/
        } else if (seq.fh[0].name) {
            EPR("== Single-end alignment");
            if ((c = map_fq_se(&seq)) != 0) {
                EPR("ERROR: map_fq_se() returned %d", c);
                goto out;
            }
            fputc('\n', stderr);
            //fq_print(&seq);
        } else {
//            if (seq.readlength == 0) {
//                EPR("== Readlength needed for indexing.");
//                goto out;
//            }
            c = fa_index(&seq);
            if (c < 0) {
                EPR("== failed to create keyct.");
                goto out;
            }
            fputc('\n', stderr);
        }
    } else {
        if (seq.fh[1].name) {
            EPR("== Paired-end assembly");
            /*if ((c = pe_fq_b2(&seq)) != 0) {
                EPR("ERROR: fq_b2() returned %d", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);*/
        } else if (seq.fh[0].name) {
            EPR("== Single-end assembly");
            if ((c = fq_b2(&seq)) != 0) {
                EPR("ERROR: fq_b2() returned %d", c);
                goto out;
            }
            fputc('\n', stderr);
            fq_print(&seq);
        } else {
            EPR("No input files specified");
            goto out;
        }
    }
    EPR("Cleanup"); fflush(NULL);

    ret = EXIT_SUCCESS;
out: /* cleanup */
    for (i=0; i != fhsz; ++i) {
        if (seq.fh[i].fp == NULL || seq.fh[i].fp == stdin || seq.fh[i].fp == stdout)
            continue;
        EPR("closing %u\n", i);
        rclose(&seq.fh[i]);
    }
    return ret;
}


