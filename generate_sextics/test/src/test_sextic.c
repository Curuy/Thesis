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

    printf("Running sextics up to 10^15...\n");

    pari_sp ltop = avma;
    
    long real_f = 0, comp_f = 0;
    // up to 10^12
    sextic_window(stoi(0), stoi(1000000000000), 0, &real_f, &comp_f);
    assert(real_f == 690);
    assert(comp_f == 2809);
    printf("[OK] 10^12 :Found %ld real, %ld complex.\n", real_f, comp_f);
    // up to 10^13
    sextic_window(stoi(0), stoi(10000000000000), 0, &real_f, &comp_f);
    assert(real_f == 1650);
    assert(comp_f == 6315);
    printf("[OK] 10^13 :Found %ld real, %ld complex.\n", real_f, comp_f);
    // up to 10^14
    sextic_window(stoi(0), stoi(100000000000000), 0, &real_f, &comp_f);
    assert(real_f == 3848);
    assert(comp_f == 14121);
    printf("[OK] 10^14 :Found %ld real, %ld complex.\n", real_f, comp_f);
    // up to 10^15
    sextic_window(stoi(0), stoi(1000000000000000), 0, &real_f, &comp_f);
    assert(real_f == 8867);
    assert(comp_f == 31276);
    printf("[OK] 10^15 :Found %ld real, %ld complex.\n", real_f, comp_f);
    avma = ltop;



    pari_close();
    return 0;
}
