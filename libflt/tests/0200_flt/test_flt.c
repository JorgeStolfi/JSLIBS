/* Quick test of rounding mode setting (from flt.h) */
/* Last edited on 2024-12-21 11:23:16 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>

#include <flt.h>

int main(int argc, char **argv);
void test_get_set_fsr(void);
void test_interp(void);
void test_int(void);

int main(int argc, char **argv)
  { flt_init();
    test_get_set_fsr();
    test_int();
    test_interp();

    fclose(stderr);
    return (0);
  }

void test_get_set_fsr(void)
  {
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, "testing ROUND_UP, DROUND_DOWN, flt_get_fsr() ...\n");

    float one, three;
    int fsr;
    one = 1.0;
    three = 3.0;

    double rdb = ((double)one)/((double)three);

    fsr = flt_get_fsr(); fprintf(stderr, "fsr   = %8x\n", fsr); 
    fprintf(stderr, "\n"); 

    ROUND_UP;              fprintf(stderr, "ROUND_UP\n"); 
    fsr = flt_get_fsr();   fprintf(stderr, "fsr   = %8x\n", fsr); 
    float rup = one/three; fprintf(stderr, "1/3   = %18.8e\n", rup);
    demand(rdb < rup, "ROUND_UP failed");
    fprintf(stderr, "\n"); 

    ROUND_DOWN;            fprintf(stderr, "ROUND_DOWN\n"); 
    fsr = flt_get_fsr();   fprintf(stderr, "fsr   = %8x\n", fsr); 
    float rdn = one/three; fprintf(stderr, "1/3   = %18.8e\n", rdn); 
    demand(rdn < rdb, "ROUND_DOWN failed");
    fprintf(stderr, "\n"); 

    double diff = ((double)rup) - ((double)rdn);   
    fprintf(stderr, "rup - rdn = %18.8e\n", diff);
    fprintf(stderr, "\n"); 

    ROUND_NEAR;            fprintf(stderr, "ROUND_NEAR\n"); 
    fsr = flt_get_fsr();   fprintf(stderr, "fsr   = %8x\n", fsr); 
    float rne = one/three; fprintf(stderr, "1/3   = %18.8e\n", rne);
    demand((rdn <= rne) && (rne <= rup), "ROUND_NEAR failed");
    fprintf(stderr, "\n"); 
  }

void test_interp(void)
  {
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, "testing flt_interp_lo, flt_interp_hi ...\n");
    /* Assumes that {Float} is {float} and not {double} !!! */
    int x0, x1, x, y0, y1;
    int dx, absy0, absy1;
    x0 = 2; 
    x1 = 5;
    dx = 3;
    absy0 = 3;
    absy1 = 1;
    for (y0 = -absy0; y0 <= +absy0; y0 += 2*absy0)
      { for (y1 = -absy1; y1 <= +absy1; y1 += 2*absy1)
          { for (x = 1; x <= 1 + 2*dx; x += dx)
              { fprintf(stderr, "interpolation between (%d,%d) and (%d,%d) at x = %d\n", x0,y0,x1,y1,x);
                /* Using {double} artithmetic: */
                double ydb = (double)y0 + ((double)x - x0)*((double)y1 - y0)/((double)x1 - x0);
                /* Using {flt_interp_lo,flt_interp_hi}: */
                Float ylo = flt_interp_lo((Float)x0,(Float)y0,(Float)x1,(Float)y1,(Float)x);
                fprintf(stderr, "lo = %24.16e  db = %24.16e  diff = %24.16e\n", ylo, ydb, ylo-ydb);
                demand(ylo <= ydb, "bug in flt_interp_lo");
                Float yhi = flt_interp_hi((Float)x0,(Float)y0,(Float)x1,(Float)y1,(Float)x);
                fprintf(stderr, "hi = %24.16e  db = %24.16e  diff = %24.16e\n", yhi, ydb, yhi-ydb);
                demand(ydb <= yhi, "bug in flt_interp_hi");
                fprintf(stderr, "\n");
              }
          }
      }
  }

void test_int(void)
  {
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, "testing flt_int ...\n");
    int i;
    srand(1234567);
    for (i = 0; i < 1000; i++)
      { int ix = (i == 0 ? INT_MIN : (i == 1 ? INT_MAX : (i == 2 ? 0 : (i <= 6 ? rand()/10000 : rand()))));
        ROUND_DOWN;
        Float fx_lo = flt_from_int(ix);
        ROUND_UP;
        Float fx_hi = flt_from_int(ix);

        bool_t bug0 = (fx_lo == -Infinity) || (fx_hi == +Infinity);
        double dx = (double)ix;
        bool_t bug1 = (((double)fx_lo > dx) || (dx > (double)fx_hi));
        double err = (double)fx_hi - (double)fx_lo;
        bool_t bug2 = (err > 1.0e-6*fabs(dx));
        bool_t bug = bug0 || bug1 || bug2;
        bool_t debug = (i < 10);
        if (bug || debug)
          { fprintf(stderr, "ix = %12d  flt_int = [ %14.7e __ %14.7e ]", ix, fx_lo, fx_hi);
            fprintf(stderr, " err = %24.16e\n", err);
            if (bug) 
              { fprintf(stderr, "** bug %c %c %c\n", "FT"[bug0], "FT"[bug1], "FT"[bug2]);
                return;
              }
          }
      }
  }
