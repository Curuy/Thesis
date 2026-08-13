#include "cubic.h"

#define DEBUG 0
#include "debug.h"

void sieve_sextic_window(GEN X, GEN delta);

void sieve_sextic_window(GEN X, GEN delta)
{
    pari_sp ltop = avma;
    // X + delta
    GEN X_delta = addii(X, delta);
    
    debug("START OF SIEVE: [%s, %s]", GENtostr(X), GENtostr(X_delta));

    // floor( ((X + delta) / 27) ^ (1/4) )
    GEN X_delta_27 = gdivent(X_delta, stoi(27));
    GEN max_f_gen = sqrtnint(X_delta_27, 4);
    long max_f = itos(max_f_gen);
    debug("Max f: %ld", max_f);
    
    long valids = 0;
    long candidates = 0;
    long real_fields = 0;
    long complex_fields = 0;    
    for (long f = 1; f <= max_f; f++)
    {
        pari_sp av_f = avma;
        GEN f_gen = stoi(f);
        GEN f4 = mulii(sqri(f_gen), sqri(f_gen));
        
        // low bound for |D2| >= (X / f^4)^(1/3)
        GEN num_min = addii(X, subis(f4, 1));
        GEN q_min = gdivent(num_min, f4);
        GEN d2_min_floor_gen = sqrtnint(q_min, 3);
        long D2_min = itos(d2_min_floor_gen);
        
        // if floor(d2^3) < q_min, incremenet
        if (cmpii(powis(d2_min_floor_gen, 3), q_min) < 0) {
            D2_min++;
        }
        
        // upper bound for |D2| <= ((X + delta) / f^4)^(1/3)
        GEN q_max = gdivent(X_delta, f4);
        GEN d2_max_floor_gen = sqrtnint(q_max, 3);
        long D2_max = itos(d2_max_floor_gen);

        // search for such a D2
        for (long d2_abs = D2_min; d2_abs <= D2_max; d2_abs++)
        {
            pari_sp av_d2 = avma; 
            
            int signs[2] = {1, -1};
            for (int s = 0; s < 2; s++)
            {
                candidates++;
                long D2 = signs[s] * d2_abs;
                GEN D2_gen = stoi(D2);
                
                if (isfundamental(D2_gen))
                {
                    GEN D3 = mulii(D2_gen, sqri(f_gen));
                    
                    if (issquareall(D3, NULL)) continue;
                    
                    valids++;
                    debug("MATCH: f = %ld, D2 = %ld, D3 = %s", f, D2, GENtostr(D3));
                    unsigned long old_HH = HH;
                    if (signe(D3) > 0) 
                    {
                        r_main(D3, D3, 0);
                        real_fields += HH - old_HH;

                    }
                    else 
                    {
                        GEN abs_D3 = negi(D3);
                        c_main(abs_D3, abs_D3, 0);
                        complex_fields += HH - old_HH;
                    }
                }
            }
            avma = av_d2;
        }
        avma = av_f;
    }
    debug("SIEVE DONE: %ld / %ld", valids, candidates);
    debug("real:%ld  complex:%ld", real_fields, complex_fields);
    avma = ltop;
}



void init_primes(unsigned long lim);
extern unsigned long HH, HH1;

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s <X> <delta>\n", argv[0]);
        exit(1);
    }
    

    pari_init(10000000, 2); 
    

    init_primes(100000); 
    HH = 0;
    HH1 = 0;
    
    GEN X = gp_read_str(argv[1]);
    GEN delta = gp_read_str(argv[2]);

    pari_timer T;
    TIMERstart(&T);

    sieve_sextic_window(X, delta);

    pari_printf("Execution Time: %ld ms\n", TIMER(&T));
    return 0;
}
