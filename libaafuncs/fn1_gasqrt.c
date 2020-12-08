/* See {fn1_gasqrt.h}. */
/* Last edited on 2016-12-26 21:27:21 by stolfilocal */
 
#include <fn1_gasqrt.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_gasqrt_eval_fp (Float x);
Interval fn1_gasqrt_eval_ia (Interval x);
AAP fn1_gasqrt_eval_aa (AAP x);
 
char *fn1_gasqrt_tag = "gasqrt";
char *fn1_gasqrt_descr = "g(x) = alt_sqrt(x)";

Float fn1_gasqrt_eval_fp (Float x)
  {
    ROUND_NEAR;
    if (x >= 0) return (Float)(sqrt(x));
    else return (Zero);
  }

Interval fn1_gasqrt_eval_ia (Interval x)
  {
    if (x.hi >= 0) return (ia_sqrt(x));
    else return (ia_const(Zero, Zero));
  }

AAP fn1_gasqrt_eval_aa (AAP x)
  {
    Interval r = aa_range(x);
    if (r.hi >= 0) return (aa_alt_sqrt(x));
    else return (aa_zero());
  }

Interval fn1_gasqrt_xd = {-One, Three};
Interval fn1_gasqrt_yd = {-Half, Two};

Float fn1_gasqrt_epsilon = 1.0e-6f;
Float fn1_gasqrt_delta = 1.0e-20f;
int fn1_gasqrt_nsub = 16;
 
fn1_data_t fn1_gasqrt_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_gasqrt_tag;
    f.descr = fn1_gasqrt_descr;
    f.eval_fp = &fn1_gasqrt_eval_fp;
    f.eval_ia = &fn1_gasqrt_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_gasqrt_eval_aa;
    f.xd = fn1_gasqrt_xd;
    f.yd = fn1_gasqrt_yd;
    f.epsilon = fn1_gasqrt_epsilon;
    f.delta = fn1_gasqrt_delta;
    f.nsub = fn1_gasqrt_nsub;
    return f;
  }
 
