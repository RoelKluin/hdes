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

int
map_fq_se(struct seqb2_t* seq)
{
    // 1) open keyindex, infior and strand
    struct gzfh_t* fhio[3] = { seq->fh + 1, seq->fh + 2, seq->fh + 3};
    const char* ext[4] = {".kc",".2b",".bd", ".fa"};
    char file[768];
    int res = -ENOMEM;
    kct_t kc = {0};
    kc.ext = seq->readlength - KEY_WIDTH;
    ASSERT(fhio[1]->name != NULL, return -EFAULT);
    unsigned len = strlen(fhio[1]->name) + 1;
    char* f = strstr(fhio[1]->name, ext[3]);
    if (f != NULL)
        strncpy(f, ext[1], strlen(ext[1]));
    else
        ASSERT(strstr(fhio[1]->name, ext[1]), return -EFAULT);

    for (int i=0; i != 3; ++i) {
        if (fhio[i]->name == NULL) {
            fhio[i]->name = &file[len*i];
            strncpy(fhio[i]->name, fhio[1]->name, len);
            _ACTION0(reopen(fhio[i], ext[0], ext[i]), 
                    "%s %s", ext[i], fhio[i]->fp ? "exists" : "does not exist")
        }
    }
    ASSERT(fhio[0]->fp && fhio[1]->fp && fhio[2]->fp, return -EFAULT,
            "need seqb2, keycount and unique boundary files");
    kc.kctndx = _buf_init_arr_err(kc.kctndx, KEYNT_BUFSZ_SHFT, return -ENOMEM);
    // 2) open seqb2 for verification of reads
    // 3) open original boundaries
    _ACTION(load_seqb2(fhio[0], &kc), "loading twobit sequence file")
    _ACTION(load_kc(fhio[1], &kc), "loading keycounts file")
    _ACTION(load_boundaries(fhio[2], &kc), "loading boundary file")
    
    // 4) open fq for reading
    //...
    EPR("All seems fine.");
err:
    EPQ(res, "an error occured:%d", res);
    free_kc(&kc);
    return res;
}
