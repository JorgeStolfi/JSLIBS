/* See {fn1_gsqrt.h}. */
/* Last edited on 2005-09-25 15:48:13 by stolfi */
 
#include <fn1_gsqrt.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_gsqrt_eval_fp (Float x);
Interval fn1_gsqrt_eval_ia (Interval x);
AAP fn1_gsqrt_eval_aa (AAP x);
 
char *fn1_gsqrt_tag = "gsqrt";
char *fn1_gsqrt_descr = "g(x) = sqrt(x)";

Float fn1_gsqrt_eval_fp (Float x)
  {
    ROUND_NEAR;
    if (x >= 0) return (sqrt(x));
    else return (Zero);
  }

Interval fn1_gsqrt_eval_ia (Interval x)
  {
    if (x.hi >= 0) return (ia_sqrt(x));
    else return (ia_const(Zero, Zero));
  }

AAP fn1_gsqrt_eval_aa (AAP x)
  {
    Interval r = aa_range(x);
    if (r.hi >= 0) return (aa_sqrt(x));
    else return (aa_zero());
  }

Interval fn1_gsqrt_xd = {-One, Three};
Interval fn1_gsqrt_yd = {-Half, Two};

Float fn1_gsqrt_epsilon = 1.0e-6;
Float fn1_gsqrt_delta = 1.0e-20;
int fn1_gsqrt_nsub = 16;
 
fn1_data_t fn1_gsqrt_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_gsqrt_tag;
    f.descr = fn1_gsqrt_descr;
    f.eval_fp = &fn1_gsqrt_eval_fp;
    f.eval_ia = &fn1_gsqrt_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_gsqrt_eval_aa;
    f.xd = fn1_gsqrt_xd;
    f.yd = fn1_gsqrt_yd;
    f.epsilon = fn1_gsqrt_epsilon;
    f.delta = fn1_gsqrt_delta;
    f.nsub = fn1_gsqrt_nsub;
    return f;
  }
 
