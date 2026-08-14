#include "cubic.h"
#include "debug.h"
#include <omp.h>

void sextic_window(GEN Y, GEN X, int num_threads, long* out_real, long* out_complex);
void init_primes(unsigned long lim);
extern __thread unsigned long HH, HH1;


