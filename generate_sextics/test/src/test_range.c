#include "sextic.h"
#include <assert.h>
#include <stdio.h>


extern __thread unsigned long HH, HH1;

extern void sextic_window(GEN X, GEN delta, int num_threads, long *out_real, long *out_complex);

int main(int argc,char **argv)
{

    pari_init(100000000, 2);
    init_primes(100000);
    HH = 0;
    HH1 = 0;

    printf("Running sextics in ranges...\n");

    pari_sp ltop = avma;
    
    long real_f = 0, comp_f = 0;
    // from 10^12 to 10^13
    sextic_window(stoi(1000000000000), stoi(10000000000000), 0, &real_f, &comp_f);
    assert(real_f == 960);
    assert(comp_f == 3506);
    printf("[OK] [10^12, 10^13] :Found %ld real, %ld complex.\n", real_f, comp_f);
    // from 10^12 to 10^14
    sextic_window(stoi(1000000000000), stoi(100000000000000), 0, &real_f, &comp_f);
    assert(real_f == 3158);
    assert(comp_f == 11312);
    printf("[OK] [10^12, 10^14] :Found %ld real, %ld complex.\n", real_f, comp_f);
    // search for specific field
    sextic_window(stoi(994011992000), stoi(994011992000), 0, &real_f, &comp_f);
    assert(real_f == 1);
    assert(comp_f == 0);
    printf("[OK] [994011992000, 994011992000] :Found %ld real, %ld complex.\n", real_f, comp_f);
    
    avma = ltop;



    pari_close();
    return 0;
}
