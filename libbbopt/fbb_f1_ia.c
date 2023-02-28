/* Defines a target function F for bbopt1.c */
/* Last edited on 2023-02-20 06:42:52 by stolfi */
  
#define _GNU_SOURCE
#include <stdint.h>

#include <ia.h>

#include <bbgoal.h>
#include <fbb_f1_ia.h>

/* PROTOTYPES */

int32_t fbb_f1_ia_dim = 1;
char *fbb_f1_ia_tag = "f1_ia";
char *fbb_f1_ia_descr = "(IA) u = x - 2/7; u^2 + 1";

Float fbb_f1_ia_eval_fp(Float *x);
Interval fbb_f1_ia_eval_ia(Interval *xr);
void fbb_f1_ia_true_opt(Interval *xr, Interval *sr);
Interval fbb_f1_ia_plot_range(Interval *xr);

/* IMPLEMENTATIONS */

#define fbb_f1_XMINN (2.0)
#define fbb_f1_XMIND (7.0)
#define fbb_f1_XMIN (fbb_f1_XMINN/fbb_f1_XMIND)

Float fbb_f1_ia_eval_fp(Float *x)
  { Float dx = (Float)((*x) - fbb_f1_XMIN);
    return (Float)(dx*dx + 1.0);
  }

Interval fbb_f1_ia_eval_ia(Interval *xr)
  { Interval x = *xr;
    Interval m = ia_scale(ia_const(One, Zero), fbb_f1_XMINN, fbb_f1_XMIND); 
    Interval dx = ia_sub(x, m);
    Interval xx = ia_sqr(dx);
    Interval f = ia_shift(xx, 1.0);
    return f;
  }

void fbb_f1_ia_true_opt(Interval *xr, Interval *sr)
  { Interval mr = ia_scale(ia_const(One, Zero), fbb_f1_XMINN, fbb_f1_XMIND); 
    Float slo, shi;
    if (mr.hi <= xr->lo)
      { slo = shi = xr->lo; }
    else if (mr.lo >= xr->hi)
      { slo = shi = xr->hi; }
    else
      { slo = (mr.lo < xr->lo ? xr->lo : mr.lo);
        shi = (mr.hi > xr->hi ? xr->hi : mr.hi);
      }
    *sr = (Interval) { slo, shi };
  }

Interval fbb_f1_ia_plot_range(Interval *xr)
  { Interval pr = fbb_f1_ia_eval_ia(xr);
    return (Interval){ 0.0, pr.hi };
  }

bbgoal_data_t fbb_f1_ia_get_data(void)
  {
    bbgoal_data_t f;
    f.dim = fbb_f1_ia_dim;
    f.eval_fp = &fbb_f1_ia_eval_fp;
    f.eval_ia = &fbb_f1_ia_eval_ia;
    f.true_opt = &fbb_f1_ia_true_opt;
    f.plot_range = &fbb_f1_ia_plot_range;
    f.tag = fbb_f1_ia_tag;
    f.descr = fbb_f1_ia_descr;
    return f;
  }

