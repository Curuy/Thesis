#if 1
#  ifndef INLINE
#    define INLINE __inline__ static
#  endif
#else /* for profiling */
#  define INLINE static
#endif

/* return (disc - Zmodif) >> shift */
INLINE ulong
get_disc_num(GEN disc)
{
  const long l = lgefint(disc)-1;
#if SHIFT
  ulong num = ((ulong)disc[l]) >> SHIFT;
  if (l > 2) num |= ((ulong)disc[l-1]) << (BITS_IN_LONG - SHIFT);
  return num - Zmodiflow;
#else
  return ((ulong)disc[l]) - Zmodiflow;
#endif
}

INLINE GEN
sqrti_up(GEN D)
{
  GEN r, sqrtD = sqrtremi(D, &r);
  return (r == gen_0)? sqrtD: addis(sqrtD, 1);
}
#define sqrti_down(D) sqrti(D)
 
/* solve quadratic equation a x^2 - 2b x + c = 0, of reduced discriminant D >= 0
 * round towards center */
INLINE void
quad_eq_c(long a, GEN b, GEN D, long *r1, long *r2)
{
  GEN sqrtD = sqrti_down(D);
  long r;
  *r1 = itos(divis_rem(subii(b,sqrtD), a, &r)); if (r > 0) (*r1)++;
  *r2 = itos(divis_rem(addii(b,sqrtD), a, &r)); if (r < 0) (*r2)--;
}
/* same, round towards infinity */
INLINE void
quad_eq_i(long a, GEN b, GEN D, long *r1, long *r2)
{
  GEN sqrtD = sqrti_up(D);
  long r;
  *r1 = itos(divis_rem(subii(b,sqrtD), a, &r)); if (r < 0) (*r1)--;
  *r2 = itos(divis_rem(addii(b,sqrtD), a, &r)); if (r > 0) (*r2)++;
}

/* modify *xmax such that
 * [xmin, xmaxnew] = [xmin, xmaxold] \cap { x,  Q(x) >= 0 }
 * return 0 if not possible, 1 otherwise
 * xmaxnew := min(xmaxold, floor((b - sqrtD) / a)) >= xmin 
 * Hence possible iff b - xmin a >= sqrtD >= 0 */
INLINE int
quad_eq1_d(long a, GEN b, GEN D, long xmin, long *xmax)
{
  GEN sqrtD, t = subii(b, mulss(a,xmin));
  long s, r, r1;
  if (signe(t) < 0 || ( s = cmpii(sqri(t), D) ) < 0) return 0;
  if (!s) { *xmax = xmin; return 1; }
  sqrtD = sqrti_up(D);
  r1 = itos(divis_rem(subii(b,sqrtD), a, &r)); if (r < 0) r1--;
  if (r1 < *xmax) *xmax = r1;
  return 1;
}
/* modify *xmin such that
 * [xminnew, xmax] = [xminold, xmax] \cap { x,  Q(x) >= 0 }
 * return 0 if not possible, 1 otherwise
 * xminnew := min(xminold, ceil((b + sqrtD) / a)) <= xmax 
 * Hence possible iff b - xmax a <= - sqrtD <= 0 */
INLINE int
quad_eq2_u(long a, GEN b, GEN D, long *xmin, long xmax)
{
  GEN sqrtD, t = subii(b, mulss(a,xmax));
  long s, r, r2;
  if (signe(t) > 0 || ( s = cmpii(sqri(t), D) ) < 0) return 0;
  if (!s) { *xmin = xmax; return 1; }
  sqrtD = sqrti_up(D);
  r2 = itos(divis_rem(addii(b,sqrtD), a, &r)); if (r > 0) r2++;
  if (r2 > *xmin) *xmin = r2;
  return 1;
}

/* binary "gcd" (strips down powers of 2 !!): x and y non-negative. */
ulong ugcd(ulong x, ulong y);
INLINE ulong
gcd(ulong a,ulong b)
{ 
  if (!b) return a;
  while (!(b&1)) b>>=1;
  return ugcd(a,b);
}

INLINE int
mod9(long x)
{
  long t = x % 9;
  if (t < 0) t += 9;
  return t;
}

INLINE void
init_U3_check(long a, long b, long c)
{
  const int bmod3 = b % 3;
  const int cmod3 = c % 3;
  _3dividesf = (bmod3 == 0 && cmod3 == 0);
  if (_3dividesf)
  {
    const int amod3 = a % 3;
    if (amod3 == 0)
    {
      forbidden_dmod9_1 = 3;
      forbidden_dmod9_2 = 6;
    }
    else
    {
      forbidden_dmod9_1 = mod9(-a-b-c);
      forbidden_dmod9_2 = mod9( a-b+c);
    }
  }
  else
  {
    if (bmod3 == 0)
      forbidden_dmod9_1 = -1; /* no restriction */
    else
    {
      const int cmod9 = c % 9;
      const int t = (bmod3 == -1 || bmod3 == 2)? 2: -2; /* 1/4b^3 mod 9 */
      forbidden_dmod9_1 = mod9(t * cmod9 * cmod9 * ((b*b - 4*a*c) % 9));
    }
  }
}

#ifdef ONLY_FIELDS
#   define INIT_U3_CHECK(a,b,c) (init_U3_check(a,b,c))
#else
#   define INIT_U3_CHECK(a,b,c)
#endif

#define _d_LOOP(f,Dmin,Dmax) { int dmod9 = mod9(Dmin); \
  DBG(if (Dmax >= Dmin) INDEX_d += Dmax - Dmin + 1;)\
  for (d=Dmin; d<=Dmax; d++) {\
    if (d && dmod9 != forbidden_dmod9_1\
          && (!_3dividesf || (dmod9 && dmod9 != forbidden_dmod9_2)))\
      f(a,b,c,d,P4,P,bc-a9*d,subis(cc,b3*d));\
    if (++dmod9 == 9) dmod9 = 0;\
  }\
}

#ifdef PARI_KERNEL_GMP
#include <gmp.h>
#define LIMBS(x)  ((mp_limb_t *)((x)+2))

/* assume x >= y */
INLINE GEN
subiuspec(GEN x, ulong s, long nx)
{
  GEN zd;
  long lz;

  if (nx == 1) return utoi(x[0] - s);

  lz = nx + 2; zd = cgeti(lz);
  mpn_sub_1 (LIMBS(zd), (mp_limb_t *)x, nx, s);
  if (! zd[lz - 1]) { --lz; }

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

/* assume x > y */
INLINE GEN
subiispec(GEN x, GEN y, long nx, long ny)
{
  GEN zd;
  long lz;
  if (ny==1) return subiuspec(x,*y,nx);
  lz = nx+2; zd = cgeti(lz);

  mpn_sub (LIMBS(zd), (mp_limb_t *)x, nx, (mp_limb_t *)y, ny);
  while (lz >= 3 && zd[lz - 1] == 0) { lz--; }

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}
#else
INLINE GEN
subiuspec(GEN x, long s, long nx)
{
  GEN xd, zd = (GEN)avma;
  long lz;
  LOCAL_OVERFLOW;

  lz = nx+2; (void)new_chunk(lz);
  xd = x + nx;
  *--zd = subll(*--xd, s);
  if (overflow)
    for(;;)
    {
      *--zd = ((ulong)*--xd) - 1;
      if (*xd) break;
    }
  if (xd == x)
    while (*zd == 0) { zd++; lz--; } /* shorten z */
  else
    do  *--zd = *--xd; while (xd > x);
  *--zd = evalsigne(1) | evallgefint(lz);
  *--zd = evaltyp(t_INT) | evallg(lz);
  avma=(long)zd; return zd;
}

INLINE GEN
subiispec(GEN x, GEN y, long nx, long ny)
{
  GEN xd,yd,zd;
  long lz;
  LOCAL_OVERFLOW;

  if (ny==1) return subiuspec(x,*y,nx);
  zd = (GEN)avma;
  lz = nx+2; (void)new_chunk(lz);
  xd = x + nx;
  yd = y + ny;
  *--zd = subll(*--xd, *--yd);
  while (yd > y) *--zd = subllx(*--xd, *--yd);
  if (overflow)
    for(;;)
    {
      *--zd = ((ulong)*--xd) - 1;
      if (*xd) break;
    }
  if (xd == x)
    while (*zd == 0) { zd++; lz--; } /* shorten z */
  else
    do  *--zd = *--xd; while (xd > x);
  *--zd = evalsigne(1) | evallgefint(lz);
  *--zd = evaltyp(t_INT) | evallg(lz);
  avma=(long)zd; return zd;
}
#endif
