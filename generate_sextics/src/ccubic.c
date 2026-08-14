/* $Id: ccubic.c,v 1.15 2001/09/28 11:22:42 kb Exp $ */
#include "cubic.h"
/*
 *  Maximum intermediate value : |l1| < 1296 (X/3)^(3/2) < 2^63 -1 
 *                                                    (with 64 bit longs)
 *        ====> X < 110,996,956,716 */

#define d_LOOP(Dmin,Dmax) _d_LOOP(c_isfield,Dmin,Dmax)

/* Z <= - disc <= X */
int
c_main(GEN Z, GEN X, long mult) 
{
  double eps = 1e-10, sq3=sqrt(3.0);
  double A,B,B0,C,C0,U,V,inva;
  long Dmin,Dmax,Cmin,Cmax;
  ulong INDEX_c = 0, INDEX_d = 0;
  long u,v,a,b,c,d,a3,a4,a9,aa27,b3,bb,bb2,bc,P4,P;
  long root2,root3,root4;

  GEN Yaa27, Zaa27;
  GEN cc,root1,T,l1,l2,l3,g1,g2;
  long sqX;

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

  if (DEBUGLEVEL) err_printf("COMPLEX fields\n");
  TAB = get_disc_array(Z,X,shift, &INC, &nbcell);
#if 0
  for (a=0; a<MAXGLOB; a++) GLOB[a] = 0;
#endif
  X = subis(Z,1);
  Z = subii(Z,INC);
  for(;;)
  {
    long ava,avb;

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
 
    B = 2*sqrt(sqX)/sqrt(sq3)+eps; A = B/sq3;
    B0 = sqX/sq3+eps;
    C0 = cbrt(gtodouble(gmul2n(X,-2)))+eps;

    avb = avma; b = bc = b3 = 0;
    for (a=1; a<A; a++) 
    {
#ifdef ONLY_FIELDS
      if (!issquarefree(stoi(a))) continue;
#endif
      a3 = 3*a; a4 = a<<2; a9 = 3*a3; aa27 = a3*a9;
      C = C0/cbrt((double) a);
      for (c=1; c<C; c++,avma=avb)
      {
        long u, r;
        INIT_U3_CHECK(a,0,c);
#ifdef ONLY_UNRAM
        if (_3dividesf || cgcd(a,c) != 1) continue;
#endif
        DBG(INDEX_c++);
        cc = mulss(c,c); P = -a3*c; P4 = P<<2;
        u = a*(a-c);
        if (u>0) {root2=sqrt((double) u)+1;} else root2=1;

        l3=mulis(cc,a4*c); l1=divis(subii(X,l3),aa27);
        root1=sqrtint(l1);

        l2=divis_rem(subii(Z,l3),aa27, &r); if (r) l2=addis(l2,1);
        if (signe(l2)>0) root2 = maxss(mysqrt(l2),root2);
        Dmax = minss(itos(root1),a+c-1);
        d_LOOP(root2,Dmax);
      }
    }
    if (DEBUGLEVEL>1) dbg_time("b",0,0,INDEX_c,INDEX_d);
 
    ava = avma;
    for (a=1; a<A; a++,avma = ava) 
    {
      long l1,root1;
#ifdef ONLY_FIELDS
      long core_a = itos((GEN)core2(stoi(a))[2]);
#endif
      a3 = 3*a; a4 = a<<2; a9 = 9*a; aa27 = a3*a9; inva = 1.0/a;
      B = a3/2.0+sqrt(B0-a3*a/4.0);
      C = C0/cbrt((double)a);
      Yaa27 = mulis(Y,aa27);
      Zaa27 = mulis(Z,aa27);

      avb = avma;
      for (b=1;b<B;b++,avma=avb)
      {
#ifdef ONLY_UNRAM
        long gcd_ab;
#endif
#ifdef ONLY_FIELDS
        if (cgcd(core_a, b) != 1) continue;
#endif
#ifdef ONLY_UNRAM
        gcd_ab = cgcd(a,b);
#endif
        bb = b*b; bb2 = bb<<1; b3 = 3*b; U = b-a; V = b+a; 

        Cmin = 1-b;
        if (b <= a3>>1)
          Cmax = C+bb*inva/3; 
        else
          Cmax = C+b-(a3/4.0);
        for (c = Cmin; c<=Cmax; c++,avma=avb) 
        {
          INIT_U3_CHECK(a,b,c);
#ifdef ONLY_UNRAM
          if (_3dividesf) continue;
          if (cgcd(gcd_ab,c) != 1) continue;
#endif
          DBG(INDEX_c++);
          P = bb-a3*c; P4 = P<<2;
          T = mulss(a9*c-bb2, b);
          g2 = addii(mulsi(P4, mulss(P,P)), Zaa27); /* 4P^3 + 27a^2 Z */
          g1 = addii(g2,Yaa27); /* 4P^3 + 27a^2 X */

          quad_eq_c(aa27,T,g1, &root1,&root2);
          Dmin = maxss( myceil((c-U)*U*inva), root1 ); 
          Dmax = minss(myfloor((c+V)*V*inva), root2 ); 
          if (Dmin > Dmax) continue;

          bc = b*c; cc = mulss(c,c); 
          if (signe(g2)>=0)
          {
            quad_eq_i(aa27,T,g2, &root1,&root2);
            root1 = minss(root1,Dmax);
            root2 = maxss(root2,Dmin);
            l1 = bb+a4*(a-c); 
            if (l1>=0)
            {
              u=mysqrt2(l1); l1=(b-u)>>1; u+=b; v=u>>1; if (u&1) v++;
              root3 = minss(l1,Dmax);
              root4 = maxss(v,Dmin);

              d_LOOP(root2,root3);
              d_LOOP(root4,root1);
              root1 = minss(root1,root3);
              root2 = maxss(root2,root4);
            } 
            d_LOOP(Dmin,root1);
            d_LOOP(root2,Dmax);
          }
          else 
          {
            l1 = bb+a4*(a-c); 
            if (l1>=0)
            {
              u = mysqrt2(l1); l1=(b-u)>>1; u+=b; v=u>>1; if (u&1) v++;
              root1 = minss(l1,Dmax);
              root2 = maxss(v,Dmin);
              d_LOOP(Dmin,root1);
              d_LOOP(root2,Dmax);
            } else d_LOOP(Dmin,Dmax);
          }
        }
      }
      if (DEBUGLEVEL>1) dbg_time("a",a,(long)A,INDEX_c,INDEX_d);
    }
    if (DEBUGLEVEL) err_printf("Time:%ld ms\n\n",timer());
#if defined(ONLY_UNRAM) && defined(COUNT_MULT)
    for(a=0; a<=nbcell; a++)
    {
      switch(TAB[a])
      {
        case 1:   GLOB[1]++; break;
        case 4:   GLOB[2]++; break;
        case 13:  GLOB[3]++; break;
        case 40:  GLOB[4]++; break;
        case 121: GLOB[5]++; break;
      }
      if (TAB[a] > 12)
      {
      /* smallest multiple of 1<<shift which is <= Z */
        ulong av = avma;
        T = T_from_shift(Z, a, shift);
	err_printf("%Ps: %d\n",negi(T),TAB[a]);
        avma = av;
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
	err_printf("%Ps: %d\n",negi(T),TAB[a]);
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
