#include "sextic.h"

void sextic_window(GEN Y, GEN X, int num_threads, long* out_real, long* out_complex)
{
    pari_sp ltop = avma;
    
    debug("START OF SIEVE: [%s, %s]", GENtostr(X), GENtostr(X_delta));

    // floor( ((X) / 27) ^ (1/4) )
    GEN X_27 = gdivent(X, stoi(27));
    GEN max_f_gen = sqrtnint(X_27, 4);
    long max_f = itos(max_f_gen);
    debug("Max f: %ld", max_f);
    
    long valids = 0;
    long candidates = 0;
    long real_fields = 0;
    long complex_fields = 0;   

    num_threads = num_threads == 0 ? omp_get_max_threads() : num_threads;
    omp_set_num_threads(num_threads);
    debug("starting sieve with %d threads...\n", num_threads);
    
    struct pari_thread *pth = malloc(num_threads * sizeof(struct pari_thread));
    
    for (int i = 1; i < num_threads; i++) {
        pari_thread_alloc(&pth[i], 10000000, NULL);
    }
    #pragma omp parallel 
    {
        int thnum = omp_get_thread_num();
        if (thnum > 0) {
            pari_thread_start(&pth[thnum]);
        }
        
        #pragma omp for schedule(dynamic, 1) reduction(+:real_fields, complex_fields, valids, candidates)
        for (long f = 1; f <= max_f; f++)
        {
            pari_sp av_f = avma;
            GEN f_gen = stoi(f);
            GEN f4 = mulii(sqri(f_gen), sqri(f_gen));
            
            // low bound for |D2| >= (Y / f^4)^(1/3)
            GEN num_min = addii(Y, subis(f4, 1));
            GEN q_min = gdivent(num_min, f4);
            GEN d2_min_floor_gen = sqrtnint(q_min, 3);
            long D2_min = itos(d2_min_floor_gen);
            
            // if floor(d2^3) < q_min, increment
            if (cmpii(powis(d2_min_floor_gen, 3), q_min) < 0) {
                D2_min++;
            }
            
            // upper bound for |D2| <= ((X) / f^4)^(1/3)
            GEN q_max = gdivent(X, f4);
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
                        
                        // keep track of # fields went through
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
        
        if (thnum > 0) {
            pari_thread_close();
        }
    }
    
    // clean threads
    for (int i = 1; i < num_threads; i++) {
        pari_thread_free(&pth[i]);
    }
    free(pth);

    printf("SIEVE DONE: %ld / %ld\n", valids, candidates);

    *out_real = real_fields;
    *out_complex = complex_fields;
    avma = ltop;
}

