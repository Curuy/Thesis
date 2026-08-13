/* $Id: rcubic.c,v 1.28 2001/09/28 11:22:42 kb Exp $ */
#include "cubic.h"
/*
 *  Maximum intermediate value :
 *       in loop: |l1| < 6.3891 X^(3/2)
 *       in cubsolve  < 62.059 X^(3/2)
 *
 * < 2^63 -1 (64-bit machine)   ====> X < ??
 *
 * < 2^31 -1 (32-bit machine)   ====> X < ?? */

void quad_eq(long a, long b, GEN c, long *r1, long *r2);

/* round towards center */
void
quad_eq(long a, long b, GEN c, long *r1, long *r2)
{
  GEN sqrtD, D = subii(mulss(b,b), mulsi((a<<2),c)); /* > 0 */
  long r;
  sqrtD = sqrti_down(D); a<<=1;
  *r1 = itos(divis_rem(subsi(-b,sqrtD), a, &r)); if (r > 0) (*r1)++;
  *r2 = itos(divis_rem(addsi(-b,sqrtD), a, &r)); if (r < 0) (*r2)--;
}

/* x^4 + a x^3 + b x^2 */
static GEN
_eval(long a, GEN b, long x)
{
  return mulii(addii(mulss(x+a, x), b), mulss(x,x));
}

long quartsolve(long a, GEN b , GEN c, long x, long y, const int inc, const int up);

/* unique root of Q(P) = P^4 + AP^3 + BP^2 = C in [x,y].
 * rounded up iff up is set. inc = 1 iff Q increasing on [x,y] */
long
quartsolve(long a, GEN b , GEN c, long x, long y, const int inc, const int up)
{
  long ltop = avma, i;
  for(;;)
  {
    i = (x+y)>>1;
    if (i == x)
    {
      if (up) { if (cmpii(_eval(a,b, x), c) != 0) i = y; }
      else    { if (cmpii(_eval(a,b, y), c) == 0) i = y; }
      avma = ltop; return i;
    }
    if (inc) { if (cmpii(_eval(a,b, i), c) > 0) y = i; else x = i; }
    else     { if (cmpii(_eval(a,b, i), c) > 0) x = i; else y = i; }
  }
}

INLINE long
sq(long x) { long _x = x; return _x * _x; }

#define d_LOOP(Dmin,Dmax) _d_LOOP(r_isfield,Dmin,Dmax)

/* Z <= disc <= X */
int
r_main(GEN Z, GEN X, long mult)
{
  const double eps = 1e-10;
  double A,B,inva9;
  long a,b,c,d,b2,b3,a3,a4,a9,aa27,bb,bb2,bc,P4,P;
  long Dmin,Dmax,Cmin,Cmax;
  ulong INDEX_c = 0, INDEX_d = 0;

  GEN Xaa27,Zaa27,Yaa27;
  GEN cc,T,g1,g2 = NULL;
  long root1,root2, sqX;
  ulong maxcubsolve;

  pari_sp ltop;
  ulong nbcell;
  GEN INC,Y, X0 = X;
#if 0
  ulong GLOB[MAXGLOB];
#endif
#ifdef SHIFT
  int shift = mult? 0: SHIFT;
#else
  int shift = 0;
#endif

  if (DEBUGLEVEL) err_printf("REAL fields\n");
  TAB = get_disc_array(Z,X,shift, &INC, &nbcell);
#if 0
  for (a=0; a<MAXGLOB; a++) GLOB[a] = 0;
#endif
  X = subis(Z,1);
  Z = subii(Z,INC);
  for(;;)
  {
    long ava,avb, Pmin,Pmax;

    DBG(INDEX_c = INDEX_d = 0);
    X = addii(X,INC); if (cmpii(X0,X) < 0) X = X0;
    Z = addii(Z,INC); init_Zmodiflow(Z,shift);
    Y = subii(X,Z); if (signe(Y) < 0) break;
    if (DEBUGLEVEL)
    {
      err_printf("m = %Ps\n",Z);
      err_printf("M = %Ps\n",X);
    }
    if (TAB) memset(TAB,0,nbcell+1);
#ifdef FULL_CHECK
    init_TI(Z, Y);
#endif

    ltop = avma;
    sqX  = itos(sqrtint(X));
    maxcubsolve = mysqrt(X);
    A = 2*sqrt((double) sqX/27);

    avb = avma; b = bc = b3 = 0;
    for (a=1;a<A;a++)
    { /* b = 0 */
#ifdef ONLY_FIELDS
      if (!issquarefree(stoi(a))) continue;
#endif
      a3 = 3*a; a4 = a<<2; a9 = 3*a3; aa27 = a3*a9;
      Xaa27 = mulis(X,aa27);
      Pmax = s_cubsolve_d(a*a9, Xaa27, maxcubsolve);
      Cmin = - (Pmax / a3);
      /* - ceil( (Z/4a)^(1/3) ) */
      Cmax = - mceil(gtodouble(gpow(gdivgs(Z,a4), ginv(stoi(3)),3)));
      if (Cmax > -a3) Cmax = -a3;
      for (c=Cmin; c <= Cmax; c++,avma=avb)
      {
        long r;
        INIT_U3_CHECK(a,0,c);
#ifdef ONLY_UNRAM
        if (_3dividesf || cgcd(a,c) != 1) continue;
#endif
        DBG(INDEX_c++);
	cc = mulss(c,c); P = -a3*c; P4 = P<<2;

	g2 = mulis(cc,a4*c);
	g1 = negi(addii(Z,g2));
        d = -itos(sqrtint(divis(g1,aa27)));
	Dmin = maxss(c/3, d);

        g1 = negi(addii(X,g2));
        if (signe(g1) > 0)
        {
          g1 = divis_rem(g1,aa27, &r); if (r) g1 = addis(g1,1);
          root1 = -mysqrt(g1);
        } else root1 = -1;
        d_LOOP(Dmin,root1);
      }
    }
    if (DEBUGLEVEL>1) dbg_time("b",0,0,INDEX_c,INDEX_d);

    ava = avma;
    for (a=1; a<A; a++, avma = ava)
    {
      long cbrt_Zaa27ov4; /* ceil( (27 a^2 Z/4)^(1/3) ) */
      long cbrt_Xaa27ov4; /* ceil( (27 a^2 X/4)^(1/3) ) */
      long l1, aa9, r1,r2;
#ifdef ONLY_FIELDS
      long core_a = itos((GEN)core2(stoi(a))[2]);
#endif
      a3 = 3*a; a4 = a<<2; a9 = 9*a; aa27 = a3*a9; inva9 = 1.0/a9;
      aa9 = a3 * a3;
      Xaa27 = mulis(X,aa27);
      Zaa27 = mulis(Z,aa27);
      Yaa27 = mulis(Y,aa27);
      cbrt_Zaa27ov4  = s_cubsolve(0,Zaa27, maxcubsolve);
      cbrt_Xaa27ov4  = s_cubsolve(0,Xaa27, maxcubsolve);

      l1 = aa27>>2; if (a&1) l1++;
      B = (a3>>1)+sqrt((double)sqX-l1);
      avb = avma;
      for (b=1; b<=B; b++, avma = avb)
      {
#ifdef ONLY_UNRAM
        long gcd_ab;
#endif
        long inf1, inf2, Q3;
        GEN Q0, Q2;
        /* disc = 0 mod 16 ? */
#if 0
        if (DEBUGLEVEL>1) dbg_time("b",b,(long)B,INDEX_c,INDEX_d);
#endif
#ifdef ONLY_FIELDS
        if (cgcd(core_a, b) != 1) continue;
#endif
#ifdef ONLY_UNRAM
        gcd_ab = cgcd(a,b);
#endif
	bb =  b*b; bb2 = bb<<1; b2 = b<<1; b3 = 3*b;

        /* 27 a^2 Z <= 4P^3 < 27 a^2 X  */
        Pmin = cbrt_Zaa27ov4;
        Pmax = cbrt_Xaa27ov4 - 1;
        Cmin = mceil(-eps + (bb - Pmax) / (double)a3);
	Cmax = mfloor(eps + (bb - Pmin) / (double)a3);
	Cmax = minss(b-a3, Cmax);
#define INIT_D_LOOP()\
          INIT_U3_CHECK(a,b,c);\
          IF_UNRAMIFIED(if (_3dividesf) continue;)\
          IF_UNRAMIFIED(if (cgcd(gcd_ab,c) != 1) continue;)\
          DBG(INDEX_c++);\
          bc = b*c; cc = mulss(c,c); P = bb-a3*c; P4 = P<<2;\
          T = mulss(a9*c-bb2, b); /* = (3ac - 2P)b */\
          Dmin = mceil((bc-P)*inva9 - eps);\
          if (a3 + b + c > 0)\
          {\
            long r;\
            Dmax = itos(divis_rem(subis(cc,P), b3, &r));\
            if (r < 0) Dmax--;\
            if (Dmin > Dmax) continue; /* false 99.999% */\
          }\
          else\
            Dmax = mfloor((bc+P)*inva9 + eps);\
          if (Z != gen_0)\
          {\
            g2 = subii(mulsi(P4, mulss(P,P)), Xaa27); /* 4P^3 - 27a^2 X */\
            g1 = addii(g2, Yaa27);                    /* 4P^3 - 27a^2 Z */\
            quad_eq_c(aa27,T,g1, &root1,&root2);\
            if (root1 > Dmin) Dmin = root1;\
            if (root2 < Dmax) Dmax = root2;\
            if (Dmin > Dmax) continue; /* false 99.9% */\
          }
/* d in first interval (if RHS(1st) = LHS(2nd) , make RHS inequality strict)*/
#define C_D_LOOPS_1(Pmin,Pmax) if ((Pmin) <= (Pmax)) {\
        Cmin = mceil(-eps + (bb - (Pmax)) / (double)a3); \
        Cmax = mfloor(eps + (bb - (Pmin)) / (double)a3); \
        Cmax = minss(b-a3, Cmax); \
        for (avb=avma, c=Cmin; c<=Cmax; c++, avma = avb)\
        {\
          INIT_D_LOOP();\
          if (Z == gen_0) g2 = subii(mulsi(P4, mulss(P,P)), Xaa27);\
          if (!quad_eq1_d(aa27,T,g2, Dmin, &Dmax)) continue;\
          if (gequal0(g2) && equalii(mulss(Dmax,aa27), T)) Dmax--;\
          d_LOOP(Dmin,Dmax);\
        }}
/* d in second interval */
#define C_D_LOOPS_2(Pmin,Pmax) if ((Pmin) <= (Pmax)) {\
        Cmin = mceil(-eps + (bb - (Pmax)) / (double)a3); \
        Cmax = mfloor(eps + (bb - (Pmin)) / (double)a3); \
        Cmax = minss(b-a3, Cmax); \
        for (avb=avma, c=Cmin; c<=Cmax; c++, avma = avb)\
        {\
          INIT_D_LOOP();\
          if (Z == gen_0) g2 = subii(mulsi(P4, mulss(P,P)), Xaa27);\
          if (!quad_eq2_u(aa27,T,g2, &Dmin, Dmax)) continue;\
          d_LOOP(Dmin,Dmax);\
        }}

	for (c=Cmin; c<=Cmax; c++, avma = avb)
	{
          INIT_D_LOOP();
          d_LOOP(Dmin,Dmax);
	}
        /* now 4P^3 - 27 a^2 X >= 0 */
        Pmin = cbrt_Xaa27ov4;
        l1 = aa9 - bb;
        Q2 = mulss(l1,l1);         /*  (b^2 - 9a^2)^2 */
        Q3 = - ((aa9 + bb) << 1);  /* - 2 (9a^2 + b^2) */
        /* inflexion points of the quartic P^4 + Q_3 P^3 + Q_2 P^2 */
        quad_eq(2, 3*(Q3>>1), Q2, &inf1,&inf2);

        /* solutions in the first d-interval */
        if (b2 > a3) goto SECOND;

        Pmax = s_cubsolve_d(sq(b2 - a3), Xaa27, maxcubsolve);
        if (Z == gen_0) { C_D_LOOPS_1(Pmin,Pmax); }
        else
        {
          long pmin = maxss(Pmin, aa9 - bb2 + 1);
          C_D_LOOPS_1(pmin, Pmax);

          Pmax = minss(pmin-1, Pmax);
          if (Pmax < Pmin) goto SECOND;

          Q0 = mulis(Zaa27, -bb);
          /* otherwise quartic ineq. never true */
          if (cmpii(_eval(Q3, Q2, inf2), Q0) <= 0)
          {
            r1 = quartsolve(Q3, Q2, Q0, inf1, inf2,0,1);
            if (r1 <= Pmax)
            {
              r2 = quartsolve(Q3, Q2, Q0, inf2,  -Q3,1,0);
              C_D_LOOPS_1(maxss(Pmin,r1), minss(Pmax,r2));
            }
          }
        }
SECOND: /* solutions in the second d-interval */
        Pmax = s_cubsolve_d(sq(b2 + a3), Xaa27, maxcubsolve);

        if (b2 > a3 && Z != gen_0) /* trivially true otherwise */
        {
          l1 = s_cubsolve(sq(b2 - a3), Zaa27, Pmax);
          if (l1 > Pmin) Pmin = l1;
        }
        l1 = aa9 - bb;
        if (l1 > Pmin) Pmin = l1;
        if (Pmax < Pmin) continue;

        Q0 = mulis(Xaa27, -bb);
        /* otherwise quartic ineq. always true */
        if (cmpii(_eval(Q3, Q2, inf2), Q0) <= 0)
        {
          r1 = quartsolve(Q3, Q2, Q0, inf1, inf2,0,0);
          if (r1 >= Pmax) { C_D_LOOPS_2(Pmin, Pmax); }
          else
          {
            r2 = quartsolve(Q3, Q2, Q0, inf2,  -Q3,1,1);
            C_D_LOOPS_2(Pmin, r1);
            C_D_LOOPS_2(r2, Pmax);
          }
        }
        else { C_D_LOOPS_2(Pmin, Pmax); }
      }
      if (DEBUGLEVEL>1) dbg_time("a",a,(long)A,INDEX_c,INDEX_d);
    }
    if (DEBUGLEVEL) dbg_cycle();
#if defined(ONLY_UNRAM) && defined(COUNT_MULT)
    for(a=0; a<=nbcell; a++)
    {
#if 0
      switch(TAB[a])
      {
	case 13:  GLOB[3]++; break;
	case 40:  GLOB[4]++; break;
	case 121: GLOB[5]++; break;
      }
#endif
      if (TAB[a] > 12)
      {
        T = T_from_shift(Z, a, shift);
	err_printf("%Ps: %d\n",T,TAB[a]); avma=ltop;
      }
    }
#endif
#if 0
    for(a=0; a<=nbcell; a++)
    {
      if (!GLOB[TAB[a]])
      {
        GLOB[TAB[a]] = 1;
        T = T_from_shift(Z, a, shift);
	err_printf("%Ps: %d\n",T,TAB[a]);
      }
    }
#endif
    pari_flush(); avma = ltop;
  }
#if 0
  for(a=1; a<=MAXGLOB; a++)
    if (GLOB[a]) err_printf("3-rank=%ld : %ld\n",a,GLOB[a]);
#endif
  return 0;
}
