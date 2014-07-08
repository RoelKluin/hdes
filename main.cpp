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
// make && zgrep -E -o "^[ACTGN]{6}" /home/roel/dev/rproject/rmap/hg19.fa.gz | head -n 100000 | valgrind ./uqct

#include <stdint.h>
#include <iostream>
#include "fq_hash.h"

using namespace std;

int main(int argc, const char* argv[])
{
    unsigned phred_offset = 33;
    if (argc > 1) {
        phred_offset = atoi(argv[1]);
        cerr << "Phred offset set to " << phred_offset << endl;
    }
    assert((phred_offset >= 33) && (phred_offset <= 126));
    
    struct fq_hash uqct(100, phred_offset);

    string nm, sq, tmp, ql;
    while ((not uqct.is_err) && getline(cin, nm) && getline(cin, sq) && getline(cin, tmp) && getline(cin, ql))
        uqct.put((const uint8_t*)sq.c_str(), (const uint8_t*)ql.c_str(), (const uint8_t*)nm.c_str());

    uqct.dump();

    if (uqct.is_err) {
        cerr << "An error occurred\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
