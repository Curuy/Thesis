/* $Id: util.c,v 1.14 2001/01/31 16:09:42 kb Exp $ */
#include "cubic.h"

//uchar *TAB;
//int _3dividesf, forbidden_dmod9_1, forbidden_dmod9_2;
//ulong Zmodiflow, HH, HH1;

//#ifdef FULL_CHECK
//static ulong *TI = NULL, *PRIMES;
//#endif
__thread uchar *TAB;
__thread int _3dividesf, forbidden_dmod9_1, forbidden_dmod9_2;
__thread ulong Zmodiflow, HH, HH1;

#ifdef FULL_CHECK
static __thread ulong *TI = NULL;
static ulong *PRIMES; /* Intentionally NOT __thread so all threads share it read-only */
#endif
__thread ulong ClusterMult;

//static ulong ClusterMult;
/* ======= Solving equations with various rounding modes ==== */

/* 1: cubic */

/* Returns an integral approximation (rounded up) to the real solution of
 * -4x^3+Ax^2+B=0. Disc must be <0 (unchecked !)
 */
long
cubsolve(long A,GEN B, GEN max)
{
  GEN i, x = stoi(A>>2), y = max, mB = negi(B);
  long ltop = avma;
  for(;;)
  {
    i = shifti(addii(x,y),-1);
    if (equalii(i,x))
    {
      if (equalii(mulii(sqri(i), subsi(A,shifti(i,2))), mB)) y = x;
      avma = ltop; return itos(y);
    }
    if (cmpii(mulii(sqri(i), subsi(A,shifti(i,2))), mB) > 0) x = i; else y = i;
  }
}

long
s_cubsolve(long A,GEN B, long max)
{
  GEN mB = negi(B);
  long ltop = avma;
  ulong x = A>>2, y = max, i;
  for(;;)
  {
    i = (x+y)>>1;
    if (i == x)
    {
      if (equalii(mulis(mulss(i,i), A - (i<<2)), mB)) y = x;
      avma = ltop; return y;
    }
    if (cmpii(mulis(mulss(i,i), A - (i<<2)), mB) > 0) x = i; else y = i;
  }
}

/* round down */
long
cubsolve_d(long A,GEN B, GEN max)
{
  GEN i, x = stoi(A>>2), y = max, mB = negi(B);
  long ltop = avma;
  for(;;)
  {
    i = shifti(addii(x,y),-1);
    if (equalii(i,x))
    {
      if (!equalii(mulii(sqri(y), subsi(A,shifti(y,2))), mB)) y = x;
      avma = ltop; return itos(y);
    }
    if (cmpii(mulii(sqri(i), subsi(A,shifti(i,2))), mB) > 0) x = i; else y = i;
  }
}

long
s_cubsolve_d(long A,GEN B, ulong max)
{
  GEN  mB = negi(B);
  long ltop = avma;
  ulong x = A>>2, y = max, i;
  for(;;)
  {
    i = (x + y) >> 1;
    if (i == x)
    {
      if (!equalii(mulis(mulss(y,y), A - (y << 2)), mB)) y = x;
      avma = ltop; return y;
    }
    if (cmpii(mulis(mulss(i,i), A - (i << 2)), mB) > 0) x = i; else y = i;
  }
}

/* 2: quadratic*/

/* returns n st n-1 < sqrt(a) <= n */
long mysqrt(GEN a) { return itos(sqrti_up(a)); }

/* returns n st n-1 <= sqrt(a) < n,  that is n = 1 + floor(sqrt(a)) */
long
mysqrt2(long a)
{
  long n = sqrt((double)a);
  if (n*n != a) n++; else n+=2;
  return n;
}

/* 3: rounding */

/* smallest integer > x */
long
myceil(double x)
{
  long res = x; /* res = trunc(x) */
  if (x >= 0 || fabs(x - res) < 1e-20) res++;
  return res;
}

/* smallest integer >= x */
long
mceil(double x)
{
  long res = x;
  if (x >= 0 && fabs(x - res) > 1e-20) res++;
  return res;
}

/* smallest integer < x */
long
myfloor(double x)
{
  long res = x;
  if (x < 0 || fabs(x - res) < 1e-20) res--;
  return res;
}

/* smallest integer <= x */
long
mfloor(double x)
{
  long res = x;
  if (x < 0 && fabs(x - res) > 1e-20) res--;
  return res;
}

/* === Modular checks ============================================== */
#ifdef FULL_CHECK
/* init PRIME[i] = i-th prime for all 5 <= p <= lim [start at i = 2] */
void
init_primes(ulong lim)
{
  long i,j, n = 0, I = sqrt((double)lim);
  char *is_prime = gpmalloc(lim+1);

  (void)memset(is_prime,1,lim+1);
  PRIMES = (ulong *) gpmalloc((lim+1)*sizeof(long));
  for (i=2; i<=I; i++)
    if (is_prime[i])
    {
      for (j=(i<<1); j<=lim; j+=i) is_prime[j] = 0;
      PRIMES[n++] = i;
    }
  for (; i<=lim; i++)
    if (is_prime[i]) PRIMES[n++] = i;
  free(is_prime);
  PRIMES = (ulong *)gprealloc((void*)PRIMES, n*sizeof(ulong));
  PRIMES[0] = evaltyp(t_VECSMALL) | evallg(n);
}

#if 0
/* sqful[i]=1 <=> p^2 | i, for some prime p>3 */
char *
init_sqfull(ulong lim)
{
  long n, p, pp, i;

  sqfull = (char *) malloc(lim+1);
  memset(sqfull, 0, lim+1);
  if (lim < 25) return;
  for(i=2;; i++) {
    p = PRIMES[i]; n = pp = p * p;
    if (pp > lim) return;
    for (n=pp; n<=lim; n+=pp) sqfull[n] = 1;
  }
  return sqfull;
}
#endif

/* assume Y <= 2^31 */
void
init_TI(GEN Z, GEN Y)
{
  long i, l = itos(Y), n = lg(PRIMES);
  ulong av = avma, j;

  if (!TI) TI = (ulong*)gpmalloc((l+1) * sizeof(ulong));

  if (DEBUGLEVEL>1) timer2();
  for (i=0; i<=l; i++) TI[i] = 1;
  if (DEBUGLEVEL>1) msgtimer("initialize TI (1)");
  for (i=2; i<n; i++,avma = av)
  {
    ulong p = PRIMES[i];
    GEN J0, r = muliu(muluu(p,p), p);

    J0 = remii(Z, r);
    if (J0 != gen_0) J0 = subii(r, J0);
    if (is_bigint(J0) || (j = itos(J0)) > l)
    {
      if (cmpii(r, Z) >= 0) break;
      continue;
    }
    if (!is_bigint(r))
    {
      ulong R = r[2];
      for ( ; j<=l; j+=R) TI[j] = 0;
    }
    else /* adding r will make j > 2^31 > l */
      TI[j] = 0;
  }
  if (DEBUGLEVEL>1) msgtimer("initialize TI (cubes)");

  for (i=2; i<n; i++,avma = av)
  {
    ulong p = PRIMES[i];
    GEN J0, r = muluu(p,p);

    J0 = remii(Z, r);
    if (J0 != gen_0) J0 = subii(r, J0);
    if (is_bigint(J0) || (j = itos(J0)) > l)
    {
      if (cmpii(r, Z) >= 0) break;
      continue;
    }
    if (!is_bigint(r))
    {
      ulong R = r[2];
      for ( ; j<=l; j+=R) TI[j] *= p;
    }
    else /* adding r will make j > 2^31 > l */
      TI[j] *= p;
  }
  if (DEBUGLEVEL>1) msgtimer("initialize TI (squares)");
}
#endif

INLINE int
lazy_check_fails(ulong A, ulong B, GEN C, GEN disc)
{
  ulong f, g = smodis(disc, 25050025); /* 5^2 ... 13^2 */
  if (g % 25 == 0 || g % 49 == 0 || g % 121== 0 || g % 169== 0) return 1;
  f = gcd(A,B);
  if (f != 1)
  {
    g = smodis(C,f);
    if (!g || gcd(f,g) != 1) return 1;
  }
  return 0;
}

#ifdef FULL_CHECK
INLINE int
full_check_fails(ulong A,ulong B,GEN C,GEN disc)
{
  const ulong f = TI[ get_disc_num(disc) ];
  return (f != 1 && (!f || A % f || B % f || smodis(C,f)));
}
#endif

/* B^2 - 4AC */
INLINE GEN
c_disc3(long A4, ulong B, GEN C)
{
  GEN a4c = mulsi(A4,C), bb;
  if (!B) return a4c;
  bb = sqru(B);
  if (signe(a4c) <= 0) return addii(bb, negi(a4c));
  return subiispec(bb+2, a4c+2, lgefint(bb)-2, lgefint(a4c)-2);
}

/* 4AC - B^2 */
INLINE GEN
r_disc3(ulong A4, ulong B, GEN C)
{
  GEN a4c = mului(A4,C), bb;
  if (!B) return a4c;
  bb = sqru(B);
  return subiispec(a4c + 2, bb+2, lgefint(a4c)-2, lgefint(bb)-2);
}

INLINE GEN
c_disc(long A4, ulong B, GEN C)  { return divis(c_disc3(A4,B,C), 3); }
INLINE GEN
r_disc(ulong A4, ulong B, GEN C) { return divis(r_disc3(A4,B,C), 3); }

#ifdef ONLY_UNRAM
#  define U_check_fails(A,B,C,disc) lazy_check_fails(A,B,C,disc)
#else
#  define U_check_fails(A,B,C,disc) full_check_fails(A,B,C,disc)
#endif

#ifdef ONLY_UNRAM
#  define U2_TEST(A4,B,C) {\
    ulong _x = (B*B - C * A4) & 0xfUL;\
    switch (_x) /* = -3disc mod 16 = disc mod 16 = 0 mod 4 */\
    {\
      case 0: case 4: return;\
    }\
}
#else
#  define U2_TEST(A4,B,C) {\
    ulong _x = (B*B - C * A4) & 0xfUL;\
    switch (_x) /* = -3disc mod 16 = disc mod 16 = 0 mod 4 */\
    {\
      case 0: return;\
      case 4: if ((A|C0) & 1) return; /* A or C odd */\
    }\
}
#endif

#ifdef CHECK_CLUSTER
#   define OUT(D) myoutput_cluster(a,b,c,d,D);
#else
# ifdef COUNT_MULT
#   define OUT(D) if (++TAB[ get_disc_num(D) ]==0) err_printf("overflow: %Ps\n",D);
# else
#   define OUT(D) { HH++; if (!HH) { HH1++; } }
# endif
#endif

#undef OUT
#if 0 /* test of Ellenberg / Venkatesh question */
long S = 0;
#define OUT(D) { S += moebius(D); err_printf("%ld ",S); }
#endif

#define OUT(D) { HH++; if (!HH) { HH1++; } }

#ifndef ONLY_FIELDS
void
r_isfield(long a,long b,long c,long d,long A4,long A,long B,GEN C) {
  GEN D = r_disc(A4,labs(B),C); OUT(D);
}
void
c_isfield(long a,long b,long c,long d,long A4,long A,long B,GEN C) {
  GEN D = c_disc(A4,labs(B),C); OUT(D);
}
#else
/* Tests boundary points and discriminant. */
void
c_isfield(long a,long b,long c,long d, long A4, long A, long B,GEN C)
{
  GEN D;

  /* test that D mod 16 != 0 or 4 (needs B even) */
  if ((B & 1) == 0)
  {
    long s = signe(C), C0;
    if (!s) C0 = 0; else C0 = (s < 0)? -mod2BIL(C): mod2BIL(C);
    U2_TEST(A4,B,C0);
  }
  /* in U_2 and U_3 */

  D = c_disc(A4,labs(B),C); /* -3disc > 0 */
  if (U_check_fails(labs(A),labs(B),C,D)) return;
#ifdef PRINT
  bigprint(a,b,c,d);
#endif
  OUT(D);
}
void
r_isfield(long a,long b,long c,long d,long A4,long A,long B,GEN C)
{
  GEN D;

  /* test that D mod 16 != 0 or 4 (needs B even) */
  if ((B & 1) == 0)
  {
    long C0 = mod2BIL(C);
    U2_TEST(A4,B,C0);
  }
  /* in U_2 and U_3 */

  { /* Test reduction for boundary point (specific to real case) */
    if (A==B) { if (b >= labs(3*a-b)) return; }
    if (!cmpis(C,A))
    {
#ifndef ONLY_UNRAM
      long tmp;
      if (a > (tmp=labs(d))) return;
      if (a==tmp) { if (b >= labs(c)) return; }
#else
      return;
#endif
    }
  } /* until here */

  D = r_disc(A4,labs(B),C);
  if (U_check_fails(A,labs(B),C,D)) return;
#ifdef PRINT
  bigprint(a,b,c,d);
#endif
  OUT(D);
}
#endif

/* ===General ====================================================== */

void error(char *msg) { err_printf(msg);exit(1);}

void
bigprint(long a, long b, long c, long d)
{
  long P = b*b-3*a*c, Q = b*c-9*a*d;
  GEN R = subis(mulss(c,c), 3*b*d);
  GEN disc = divis(subii(mulsi(P<<2,R), mulss(Q,Q)), 3);
#ifdef TEX
  GEN form = mkpoln(4, stoi(d),stoi(c),stoi(b),stoi(a));
  pari_printf("$%Ps$ & %Ps\n",form,disc);
#else
  pari_printf("%Ps: (%ld,%ld,%Ps)  [%ld,%ld,%ld,%ld]\n",disc,P,Q,R,a,b,c,d);
#endif
  pari_flush();
}

void
myoutput_cluster(long a,long b,long c,long d,GEN disc)
{
#ifdef ONLY_UNRAM
  static int MaxReached = 0;
#endif
  ulong num;

  if (!direct_U_check(a,b,c,d)) { ClusterMult--; return; }
  num = get_disc_num(disc);
  if (++TAB[num] == 0) bigprint(a,b,c,d); /* overflow */
#ifdef ONLY_UNRAM
  switch(TAB[num])
  {
    case 1: ClusterMult--;
      if (MaxReached < 1) MaxReached = 1;
      break;
    case 2: /* >= 4  */ ClusterMult -= 3;
      if (MaxReached < 4) MaxReached = 4;
      break;
    case 5: /* >= 13 */ ClusterMult -= 9;
      if (MaxReached < 13) MaxReached = 13;
      break;
    case 14:/* >= 40 */ ClusterMult -= 27; bigprint(a,b,c,d);
      break;
  }
#if 0
  bigprint(a,b,c,d);
  err_printf("TAB[%ld] = %ld, ClusterMult = %ld, MaxReached = %ld\n",
              num,TAB[num],ClusterMult,MaxReached);
#endif
  if (ClusterMult + MaxReached < 40)
  {
    err_printf("Global Time: %ld ms\n",timer2());
    exit(0);
  }
#endif
}

/* Will look for disc, Z <= disc <= X. Return multiplicity array TAB, for a
 * subrange of INC discriminant (will need about (X-Z) / INC passes).
 * nbcell = number of elements of TAB */
uchar *
get_disc_array(GEN Z, GEN X, int shift, GEN *INC, ulong *nbcell)
{
  GEN t;
#ifdef FULL_CHECK
  *nbcell = 1 << 26; /* <--> 64*4 MB */
#else
  *nbcell = 1 << 27; /* 128MB, was 1<<25 = 32 MB */
#endif
  *INC = addis(subii(X,Z), 1);
  t = addis(shifti(*INC, -shift), 1); /* # of needed cells */
  if (cmpis(t,*nbcell) < 0)
    *nbcell = itos(t);
  else
    *INC = shifti(stoi(*nbcell), shift);
#ifdef COUNT_MULT
  return gpmalloc(*nbcell+1);
#else
  return NULL;
#endif
}

/* Let Z = largest integer Zmodif <= Z, Zmodif = 0 (2^shift).
 * Sets Zmodiflow = lowest word of Z >> shift */
void
init_Zmodiflow(GEN Z, int shift)
{
  long l = lgefint(Z)-1;
  if (l == 1) { Zmodiflow = 0; return; }
  Zmodiflow = (ulong)Z[l];
  if (shift && l > 2)
  {
    Zmodiflow >>= shift;
    Zmodiflow |= ((ulong)Z[l-1]) << (BITS_IN_LONG - shift);
  }
}
GEN
T_from_shift(GEN Z, long a, int shift)
{ /* smallest multiple of 1<<shift which is <= Z */
  GEN Zmodif = icopy(Z);
  Zmodif[lgefint(Z)-1] = Zmodiflow << shift;
  return addii(Zmodif, shifti(stoi(a),shift));
}

/* debugging output */

int
direct_U_check(long a,long b,long c,long d)
{
  long A = b*b - 3*a*c;
  long B = b*c - 9*a*d;
  GEN C = subis(mulss(c,c), 3*b*d);
  long fH, t = gcd(labs(A),labs(B));
  GEN T,D,disc;
    
  fH = gcd(t, labs(smodis(C,t)));
  disc = subii(mulss(B,B), mulsi(4*A,C));
  T = divii(disc, muluu(fH,fH));
  D = divis(disc, -3);

  t = mod16(D);
  if (!t) return 0;

  if (signe(D) < 0) t = 16 - t;
  if (t == 4 && ((A|B|mod2BIL(C))& 1)) return 0;
  /* in U2 */

  if (b % 3 == 0 && c % 3 == 0)
  {
    long a3 = a % 3, d3 = d % 3;
    if (a3 < 0) a3 += 3;
    if (d3 < 0) d3 += 3;
    if      (a3 == 0)  { if (a % 9 == 0 || d % 3 == 0) return 0; }
    else if (d3 == 0)  { if (d % 9 == 0) return 0; }
    else if (a3 == d3) { if ((a - b + c - d) % 9 == 0) return 0; }
    else                 if ((a + b + c + d) % 9 == 0) return 0;
  }
  else
    if (smodis(D,9) == 0) return 0; /* 3 \nmid fH */
  /* in U_3 */

  t = fH; while (t % 3 == 0)  t /= 3;
  if (!issquarefree(stoi(t))) return 0;
  T = shifti(T, -vali(T));
  while (smodis(T,3) == 0) T = divis(T,3);
  if (!gequal1(ggcd(T,stoi(t)))) return 0;
  return issquarefree(T);
}

void 
dbg_H()
{
#ifdef FULL_CHECK
  if (HH1)
    err_printf("H = %Ps", mkintn(2,HH1,HH));
  else
    err_printf("H = %lu", HH);
#endif
}

void
dbg_time(char *s, long a, long A, ulong ic, ulong id)
{
#ifdef FULL_CHECK
  dbg_H(); err_printf(",\t");
#endif
  msgtimer(" %s=%ld\t/%ld\t(%lu,%lu)",s,a,(long)A,ic,id);
}

void
dbg_cycle()
{
#ifdef FULL_CHECK
  dbg_H(); err_printf(",\t");
#endif
  err_printf("Time:%ld ms\n\n",timer());
}

#ifndef EXCLUDE_MAIN
int
main(int argc,char **argv)
{
  char *s;
  long i, sx, sy, CMult = 0;
  GEN X,Z;
  pari_timer T;

  pari_init(5000000,2);
  X = Z = NULL;

  for (i=1; i<argc; i++)
  {
    if (!strncmp(argv[i],"-c",2))
    { /* check cluster, [start, mult, size] */
#ifndef CHECK_CLUSTER
      error("recompile with CHECK_CLUSTER enabled");
#endif
      s = argv[i]+2;
      if (!*s)
      {
	i++;
	if (i>=argc) error("I waited for a disc");
	s=argv[i];
      }
      X = gp_read_str(s);
      if (typ(X) != t_VEC || lg(X) != 4) error("incorrect cluster input");
      Z = (GEN)X[1];
      CMult = itos((GEN)X[2]);
      X = addis(Z, itos((GEN)X[3]));
    }
    else if (!strncmp(argv[i],"-g",2)) /* debuglevel */
    {
      s = argv[i]+2;
      if (!*s)
      {
	i++;
	if (i>=argc) { DEBUGLEVEL=1; break; }
	s=argv[i];
      }
      DEBUGLEVEL = atol(s);
    }
    else X = gp_read_str(argv[i]);
  }
  if (!X) error("Need discriminant bounds\n");
  if (typ(X) == t_INT) { sx = signe(X); X = absi(X); }
  else
  {
    long lx = lg(X);
    if (typ(X) != t_VEC || lx == 1 || lx > 3) error("incorrect parameters\n");
    Z = gel(X,1); sy = signe(Z); Z = absi(Z);
    if (lx == 2)
    {
      X = Z; sx = sy;
    }
    else /* lx == 3 */
    {
      GEN Y;
      X = gel(X,2); sx = signe(X); X = absi(X);
      if (sx * sy < 0) error("need same signs");
      Y = gmax(X,Z);
      Z = gmin(X,Z); X = Y;
    }
  }
  if (!Z)  Z = gen_0;

  TIMERstart(&T);
#ifdef FULL_CHECK
  init_primes(mysqrt(X));
  HH = HH1 = 0;
#endif
  ClusterMult = CMult;

  timer(); timer2();
  if (sx < 0) c_main(Z,X,CMult); else r_main(Z,X,CMult);
  dbg_H(); err_printf("\nGlobal Time: %ld ms\n", TIMER(&T));
  return 0;
}
#endif
