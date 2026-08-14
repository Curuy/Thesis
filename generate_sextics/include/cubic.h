#include "config.h"
#include "pari.h"
typedef unsigned char uchar;

#if !defined(COUNT_MULT) || defined(CHECK_CLUSTER)
#  undef SHIFT
#  define SHIFT 0
#endif

#ifdef FULL_CHECK
#  define ONLY_FIELDS
#endif

#ifdef ONLY_UNRAM
#  define ONLY_FIELDS
#  define IF_UNRAMIFIED(x) x
#else
#  ifdef ONLY_FIELDS
#    define FULL_CHECK
#  endif
#  define IF_UNRAMIFIED(x)
#endif

/*   from util.c  */
int f_test_fails(long f, GEN disc3);
void dbg_H(void);
void dbg_cycle(void);
void dbg_time(char *s, long a, long A, ulong ic, ulong id);
void myoutput_cluster(long a,long b,long c,long d,GEN D);
int direct_U_check(long a,long b,long c,long d);
uchar *get_disc_array(GEN Z, GEN X, int shift, GEN *INC, ulong *nbcell);
GEN T_from_shift(GEN Z, long a, int shift);
void init_Zmodiflow(GEN Z, int shift);
char find_in_list(int *t,int x,int a,int b);
double myabs(double x);
long cubsolve(long A,GEN B, GEN max);
long cubsolve_d(long A,GEN B, GEN max);
long s_cubsolve(long A,GEN B, long max);
long s_cubsolve_d(long A,GEN B, ulong max);
long mysqrt(GEN a);
long mysqrt2(long a);
long find_P(double x,int a,int b);
long init_all(int argc,char **argv);
long mceil(double x);
long mfloor(double x);
long myceil(double x);
long myfloor(double x);
void bigprint(long a, long b, long c, long d);
void error(char *msg);
void init_TI(GEN Z, GEN Y);
int testf(long f,int a,int b,int c,int d, GEN disc3);
int testf_cyclic(long f, long a, long b, long d);

/*  from rcubic.c / ccubic.c */
int main(int argc,char **argv);
int r_main(GEN Z, GEN X, long mult);
int c_main(GEN Z, GEN X, long mult);
void r_isfield(long a,long b,long c,long d,long A4,long A,long B,GEN C);
void c_isfield(long a,long b,long c,long d,long A4,long A,long B,GEN C);
void c_myoutput(int a,int b,int c,int d,GEN disc);
void r_myoutput(int a,int b,int c,int d,GEN disc);
void r_myoutput_cluster(int a,int b,int c,int d,GEN disc);

#if 1
#  define DBG(x) x
#else
#  define DBG(x)
#endif

//extern uchar *TAB;
//extern int _3dividesf, forbidden_dmod9_1, forbidden_dmod9_2;
//extern ulong Zmodiflow, HH, HH1;
extern __thread uchar *TAB;
extern __thread int _3dividesf, forbidden_dmod9_1, forbidden_dmod9_2;
extern __thread ulong Zmodiflow, HH, HH1;


extern GEN muluu(ulong x, ulong y);
extern GEN mului(ulong x, GEN y);

#define MAXGLOB 100

#include "cubicinl.h"
