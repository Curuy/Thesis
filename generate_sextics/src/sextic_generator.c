#include "sextic.h"


int main(int argc, char** argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s <Y> <X> <n threads>\n", argv[0]);
        exit(1);
    }
    int num_threads = 0;
    pari_init(10000000, 2); 
    if (argc == 4)
    {
        // thread number specified
        GEN g = gp_read_str(argv[3]);
        if (typ(g) != t_INT) 
        {
            fprintf(stderr, "error: thread number must be integer.\n");
            pari_close();
            return 1;
        }
        num_threads = itos(g);
    }


    init_primes(100000); 
    HH = 0;
    HH1 = 0;
    
    GEN Y = gp_read_str(argv[1]);
    GEN X = gp_read_str(argv[2]);
    if (gcmp(Y, X) > 0)
    {
        fprintf(stderr, "error: Y <= X\n");
        pari_close();
        return 1;
    }
    pari_timer T;
    TIMERstart(&T);
    long real_fields = 0;
    long complex_fields = 0;
    sextic_window(Y, X, num_threads, &real_fields, &complex_fields);
    printf("real:%ld  complex:%ld, total: %ld\n", real_fields, complex_fields, real_fields + complex_fields);
    pari_printf("Execution Time: %ld ms\n", TIMER(&T));
    pari_close();
    return 0;
}
