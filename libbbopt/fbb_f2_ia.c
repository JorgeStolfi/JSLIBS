/* Defines a target function F for bbopt1.c */
/* Last edited on 2017-01-02 13:41:59 by jstolfi */
  
#include <fbb_f2_ia.h>
#include <bbgoal.h>

#include <ia.h>
#include <bool.h>
#include <affirm.h>
#include <limits.h>

/* PROTOTYPES */

int fbb_f2_ia_dim = 1;
char *fbb_f2_ia_tag = "f2_ia";
char *fbb_f2_ia_descr = "(IA) u = x - 1/7; v = x - 6/7; u^2(v^2 + 1/32) + 1/32";

Float fbb_f2_ia_eval_fp(Float *x);
Interval fbb_f2_ia_eval_ia(Interval *xr);
void fbb_f2_ia_true_opt(Interval *xr, Interval *sr);
Interval fbb_f2_ia_plot_range(Interval *xr);

/* IMPLEMENTATIONS */

#define ia_One (Interval){1.0, 1.0}

#define fbb_f2_XUN (6.0)
#define fbb_f2_XUD (7.0)

#define fbb_f2_XVN (1.0)
#define fbb_f2_XVD (7.0)

#define fbb_f2_EPS (1.0/512.0)

Float fbb_f2_ia_eval_fp(Float *x)
  { Float u = (Float)((*x) - fbb_f2_XUN/fbb_f2_XUD);
    Float v = (Float)((*x) - fbb_f2_XVN/fbb_f2_XVD);
    return (Float)(u*u*(v*v + fbb_f2_EPS) + fbb_f2_EPS);
  }

Interval fbb_f2_ia_eval_ia(Interval *xr)
  { Interval x = *xr;
    Interval mu = ia_scale(ia_One, fbb_f2_XUN, fbb_f2_XUD); 
    Interval u = ia_sub(x, mu);
    Interval uu = ia_sqr(u);
    Interval mv = ia_scale(ia_One, fbb_f2_XVN, fbb_f2_XVD); 
    Interval v = ia_sub(x, mv);
    Interval vv = ia_sqr(v);
    Interval a = ia_shift(vv, fbb_f2_EPS);
    Interval b = ia_mul(uu, a);
    Interval f = ia_shift(b, fbb_f2_EPS);
    return f;
  }
  
/* 
  Function: 
    u^2(v^2 + EPS) + EPS
  Derivative:  
    2u(v^2 + EPS) + u^2(2 v) = 
    2u(v^2 + v u + EPS) = 
    2u(v^2 + v(v + v0 - u0)) + EPS) = 
    2u(2v^2 + v D + EPS)
  where u = x - u0, v = x - v0, D = v0 - u0. 
  Assuming v0 < u0:
  First minimum at    x = v0 + (sqrt(D^2 - 8 EPS) - D)/2
  Internal maximum at x = v0 - (sqrt(D^2 - 8 EPS) + D)/2
  Second minimum at   x = u0 */

void fbb_f2_ia_true_opt(Interval *xr, Interval *sr)
  { Interval C = ia_const(fbb_f2_EPS, Zero); 
    Interval u0 = ia_scale(ia_One, fbb_f2_XUN, fbb_f2_XUD);
    Interval v0 = ia_scale(ia_One, fbb_f2_XVN, fbb_f2_XVD);
    Interval D = ia_sub(v0, u0);
    Interval sd = ia_sqrt(ia_sub(ia_sqr(D), ia_scale(C, 8.0, 1.0)));
    Interval xmin1 = ia_add(v0, ia_scale(ia_sub(sd, D), 1.0, 2.0)); 
    Interval xmax = ia_sub(v0, ia_scale(ia_add(sd, D), 1.0, 2.0));
    Interval xmin2 = u0;
    
    auto void do_sol(Interval *xr, Interval *sr);
    void do_sol(Interval *xr, Interval *sr)
      { fprintf(stderr, "  xr  = "); ia_print(stderr, *xr); fprintf(stderr, "\n");
        Float slo, shi;
        if (xr->lo > xr->hi)
          { slo = xr->lo; shi = xr->hi; }
        if (xr->lo >= xmin2.hi )
          { slo = shi = xr->lo; }
        if (xr->hi <= xmin1.lo )
          { slo = shi = xr->hi; }
        else if (xr->lo >= xmax.hi)
          { if (xr->hi >= xmin2.hi)
              { slo = xmin2.lo; shi = xmin2.hi; }
            else if (xr->hi >= xmin2.lo)
              { slo = xmin2.lo; shi = xr->hi; }
            else
              { slo = shi = xr->hi; }
          }
        else if (xr->hi <= xmax.lo)
          { if (xr->lo <= xmin1.lo)
              { slo = xmin1.lo; shi = xmin1.hi; }
            else if (xr->lo <= xmin1.hi)
              { slo = xr->lo; shi = xmin1.hi; }
            else
              { slo = shi = xr->lo; }
          }
        else
          { affirm(FALSE, "bad do_sol call"); /* GCC pacifier: */ slo = shi = NAN; }
        *sr = (Interval) { slo, shi };
      }
      
    fprintf(stderr, "u0    = "); ia_print(stderr, u0); fprintf(stderr, "\n");
    fprintf(stderr, "v0    = "); ia_print(stderr, v0); fprintf(stderr, "\n");
    fprintf(stderr, "sd    = "); ia_print(stderr, sd); fprintf(stderr, "\n");
    fprintf(stderr, "xmin1 = "); ia_print(stderr, xmin1); fprintf(stderr, "\n");
    fprintf(stderr, "xmax  = "); ia_print(stderr, xmax); fprintf(stderr, "\n");
    fprintf(stderr, "xmin2 = "); ia_print(stderr, xmin2); fprintf(stderr, "\n");
    if ((xr->lo >= xmax.hi) || (xr->hi <= xmax.lo))
      { do_sol(xr, sr); }
    else
      { /* The interval {xr} may span the maximum and hence inresect both basins. */
        /* Trisect {xr} at {xmax} and find an eclosure for the minimum in each part: */
        Interval sra, srm, srb, yra, yrm, yrb;
        Float yminhi = +INFINITY; /* Upper bound for minimum in {xr}. */
        Interval xra = ia_meet((Interval){xr->lo, xmax.lo}, *xr);
        if (xra.lo <= xra.hi) 
          { do_sol(&xra, &sra); yra = fbb_f2_ia_eval_ia(&sra); 
            if (yra.hi <= yminhi) { yminhi = yra.hi; } 
          }
        else
          { yra = (Interval){ 1.0, -1.0 }; }
        Interval xrm = ia_meet(xmax, *xr);
        if (xrm.lo <= xrm.hi) 
          { srm = xrm; yrm = fbb_f2_ia_eval_ia(&srm); 
            if (yrm.hi <= yminhi) { yminhi = yrm.hi; }
          }
        else
          { yrm = (Interval){ 1.0, -1.0 };
          }
        Interval xrb = ia_meet((Interval){xmax.hi, xr->hi}, *xr);
        if (xrb.lo <= xrb.hi) 
          { do_sol(&xrb, &srb); yrb = fbb_f2_ia_eval_ia(&srb); 
            if (yrb.hi <= yminhi) { yminhi = yrb.hi; }
          }
        else
          { yrb = (Interval){ 1.0, -1.0 }; }

        /* Join the intervals {sra, srm, srb} may contain minimum: */
        *sr = (Interval){ 1.0, -1.0 };
        if ((yra.hi >= yra.lo) && (yra.lo <= yminhi)) { *sr = ia_join(*sr, sra); }
        if ((yrm.hi >= yrm.lo) && (yrm.lo <= yminhi)) { *sr = ia_join(*sr, srm); }
        if ((yrb.hi >= yrb.lo) && (yrb.lo <= yminhi)) { *sr = ia_join(*sr, srb); }
      }
  }

Interval fbb_f2_ia_plot_range(Interval *xr)
  { Interval w0 = ia_scale(ia_One, fbb_f2_XUN + fbb_f2_XVN, fbb_f2_XUD + fbb_f2_XVD);
    Interval w = ia_sub(*xr, w0);
    return ia_shift(ia_sqr(ia_sqr(w)), fbb_f2_EPS);
  }

bbgoal_data_t fbb_f2_ia_get_data(void)
  {
    bbgoal_data_t f;
    f.dim = fbb_f2_ia_dim;
    f.eval_fp = &fbb_f2_ia_eval_fp;
    f.eval_ia = &fbb_f2_ia_eval_ia;
    f.true_opt = &fbb_f2_ia_true_opt;
    f.plot_range = &fbb_f2_ia_plot_range;
    f.tag = fbb_f2_ia_tag;
    f.descr = fbb_f2_ia_descr;
    return f;
  }

