/* See {fn1_gsqr.h}. */
/* Last edited on 2016-12-26 21:34:09 by stolfilocal */
 
#include <fn1_gsqr.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_gsqr_eval_fp (Float x);
Interval fn1_gsqr_eval_ia (Interval x);
AAP fn1_gsqr_eval_aa (AAP x);
 
char *fn1_gsqr_tag = "gsqr";
char *fn1_gsqr_descr = "g(x) = sqr(x)";

Float fn1_gsqr_eval_fp (Float x)
  {
    ROUND_NEAR;
    return(x*x);
  }

Interval fn1_gsqr_eval_ia (Interval x)
  {
    return (ia_sqr(x));
  }

AAP fn1_gsqr_eval_aa (AAP x)
  {
    return (aa_sqr(x));
  }

Interval fn1_gsqr_xd = {-Two, Two};
Interval fn1_gsqr_yd = {-One, One+Four};

Float fn1_gsqr_epsilon = 1.0e-6f;
Float fn1_gsqr_delta = 1.0e-20f;
int fn1_gsqr_nsub = 16;
 
fn1_data_t fn1_gsqr_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_gsqr_tag;
    f.descr = fn1_gsqr_descr;
    f.eval_fp = &fn1_gsqr_eval_fp;
    f.eval_ia = &fn1_gsqr_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_gsqr_eval_aa;
    f.xd = fn1_gsqr_xd;
    f.yd = fn1_gsqr_yd;
    f.epsilon = fn1_gsqr_epsilon;
    f.delta = fn1_gsqr_delta;
    f.nsub = fn1_gsqr_nsub;
    return f;
  }
 
