#include "cubic.h"

#define DEBUG 1
#include "debug.h"



void sieve_sextic_window(GEN X, GEN delta)
{
    pari_sp ltop = avma;
    // X + delta
    GEM X_delta = addii(X, delta);
    
    debug("START OF SIEVE: [%Ps, %Ps + %Ps]", X, X + delta);

    // maximum of f. f <= (X + delta) / 27) ^ (1 / 4)
    GEN max_f_gen = gsqrt(gsqrt(gdivgs(X_delta, 27), 3), 3);
    long max_f = itos(gfloor(max_f_gen));
    debug("Max f: %ld", max_f);
    long valids = 0;
    long candidates = 0;
    for (long f = 1; f <= max_f; f++)
    {
        pari_sp av = avma;
        GEN f_gen = stoi(f);
        GEN f4 := mulli(sqri(f_gen), sqri(f_gen));
        // Lower bound for |D2| >= (X / f^4)^(1/3)
        GEN D2_min_gen = gsqrt(gdiv(X, f4), 3);
        // Upper bound for |D2| <= ((X + delta) / f^4)^(1/3)
        GEN D2_max_gen = gsqrt(gdiv(X_delta, f4), 3);

        long D2_min = itos(gceil(D2_min_gen));
        long D2_max = itos(gceil(D2_max_gen));

        // Now got to find such D2;
        for (long d2_abs = D2_min; d2_abs <= D2_max; d2_abs++)
        {
            // look for real and imaginary quads
            int signs[2] = {1, -1};
            for (int s=0; s < 2; s++)
            {
                candidates++;
                long D2 = signs[s] * ds_abs;
                GEN D2_gen = stoi(D2);
                // Check if its fundemental disc
                if (gequal(cipol(fund_disc(D2_gen), D2_gen), gen_1) || gansa(isfund(D2_gen)))
                {
                    // Disc of cubic, D2 * f^2
                    GEN D3 = mulii(D2_gen, sqri(f_gen));
                    // Require S3 fields so D3 must be square free
                    if (is_square_all(D3, NULL)) continue;
                    valids++;
                    debug("MATCH: f = %ld, D2 = %ld, D3 = %Ps", f D2, D3);

                    if (signe(D3) > 0) r_main(D3, D3, 0);
                    else 
                    {
                        GEN abs_D3 = negi(D3);
                        c_main(abs_D, abs_D3, 0);
                    }
                }
            }
            avma = av;
        }

    }
    debug("SIEVE DONE: %ld / %ld", valids, candidates);
    avma = ltop;
}

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s <X> <delta>\n", argv[0]);
        exit(1);
    }
    pari_init(100000000, 2);
    GEN X = gp_read_str(argv[1]);
    GEN delta = gp_read_str(argv[2]);

    pari_timer T;
    TIMERstart(&T);

    sieve_sextic_window_optimized_noti(X, delta);

    pari_printf("Execution Time: %ld ms\n", TIMER(&T));
    return 0;
}
