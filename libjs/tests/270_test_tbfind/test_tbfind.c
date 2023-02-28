#define PROG_NAME "test_tbfind"
#define PROG_DESC "tests the ordered table search procedure"
#define PROG_VERS "1.1"

/* Last edited on 2023-02-25 16:08:01 by stolfi */
/* Created on 2003-09-25 or earler by J. Stolfi, UNICAMP */

#define test_tbfind_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

#include <tbfind.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>
#include <rn.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

int main (int argc, char **argv);
double *make_rnd_table(int n, double fMin, double fMax);
double *make_lin_table(int n, double fMin, double fMax);
double *make_exp_table(int n, double fMin, double fMax);
void test_table(char *kind, int n, double *v, int nSteps, double fMin, double fMax);

double *make_lin_table(int n, double fMin, double fMax)
  { double *v = rn_alloc(n);
    int i;
    for (i = 0; i < n; i++)
      { double r = ((double)i)/((double)(n-1));
        v[i] = (1-r)*fMin + r*fMax;
      }
    return v;
  }

double *make_rnd_table(int n, double fMin, double fMax)
  { double *v = rn_alloc(n);
    int i;
    for (i = 0; i < n; i++)
      { double s = drandom();
        double t = ((double)i)/((double)(n-1));
        double r = 0.1*s + 0.9*t;
        v[i] = (1-r)*fMin + r*fMax;
      }
    return v;
  }

double *make_exp_table(int n, double fMin, double fMax)
  { double *v = rn_alloc(n);
    double fBas = (fMin >= 0 ? 0 : fMin) - 0.001*(fMax - fMin);
    double h = log((fMax - fBas)/(fMin - fBas));
    int i;
    for (i = 0; i < n; i++)
      { double r = ((double)i)/((double)(n-1));
        v[i] = fBas + (fMin - fBas) * exp(r*h);
      }

    return v;
  }

int main (int argc, char **argv)
  {
    double fMin = -0.001, fMax = +3.000;
    int nSteps = 1000;  /* Number of levels to search: */
    srandom(4615);

    fprintf(stderr, "%s %6s", "tbl", "n");
    fprintf(stderr, " %5s %7s  %12s %7s %6s\n", "finds", "probes", "probes/fnd", "avg", "bpp");
    
    fprintf(stderr, "%s %6s", "---", "------");
    fprintf(stderr, " %5s %7s  %12s %7s %6s\n", "-----", "-------", "------------", "-------", "------");
    
    
    
    int n;
    for (n = 100; n <= 10000; n *= 10)
      { 
        double *vLin = make_lin_table(n, fMin, fMax);
        test_table("lin", n, vLin, nSteps, fMin, fMax);

        double *vExp = make_exp_table(n, fMin, fMax);
        test_table("exp", n, vExp, nSteps, fMin, fMax);

        double *vRnd = make_rnd_table(n, fMin, fMax);
        test_table("rnd", n, vRnd, nSteps, fMin, fMax);
      }
    return 0;
  }
  
void test_table(char *kind, int n, double *v, int nSteps, double fMin, double fMax)
  { int i;
    int iMin = 418, iMax = iMin + n - 1;
    int nfinds = 0, nprobes = 0;
    int npMin = INT_MAX, npMax = 0;
  
    fprintf(stderr, "%s %6d", kind, iMax-iMin+1);
    for (i = 0; i <= nSteps; i++)
      { double r = ((double)i-2)/((double)nSteps - 4); 
        double z = (1-r)*fMin + r*fMax;

        auto double f(int i);
        double f(int i)
          { nprobes++;
            affirm(i >= iMin, "i < iMin");
            affirm(i <= iMax, "i > iMax");
            return v[i - iMin] - z; 
          }

        int npOld = nprobes;
        int k =  tb_find(f, iMin, iMax); nfinds++;
        affirm((k == iMin) || (v[k-1-iMin] - z < 0), "bad k (below)");
        affirm((k == iMax+1) || (0 <= v[k-iMin] - z), "bad k (above)");
        int np = nprobes - npOld;
        if (np < npMin) { npMin = np; }
        if (np > npMax) { npMax = np; }
      }
    double ppf = ((double)nprobes)/((double)nfinds); /* Avg. probes per find. */
    double bpp = (log(n+1)/log(2))/ppf; /* Avg. bits gained pr probe. */
    fprintf
      ( stderr, " %5d %7d  %4d .. %4d %7.3f %6.3f\n", 
        nfinds, nprobes, npMin, npMax, ppf, bpp
      );
  }
  
