/* See {fn1_ginv.h}. */
/* Last edited on 2005-09-25 15:47:45 by stolfi */
 
#include <fn1_ginv.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_ginv_eval_fp (Float x);
Interval fn1_ginv_eval_ia (Interval x);
AAP fn1_ginv_eval_aa (AAP x);
 
char *fn1_ginv_tag = "ginv";
char *fn1_ginv_descr = "g(x) = 1/x";

Float fn1_ginv_eval_fp (Float x)
  {
    ROUND_NEAR;
    if (x != Zero) return (One/x);
    else return (Zero);
  }

Interval fn1_ginv_eval_ia (Interval x)
  {
    return (ia_inv(x));
  }

AAP fn1_ginv_eval_aa (AAP x)
  {
    return (aa_inv(x));
  }

Interval fn1_ginv_xd = {-Four, Four};
Interval fn1_ginv_yd = {-Four, Four};

Float fn1_ginv_epsilon = 1.0e-6;
Float fn1_ginv_delta = 1.0e-20;
int fn1_ginv_nsub = 31;
 
fn1_data_t fn1_ginv_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_ginv_tag;
    f.descr = fn1_ginv_descr;
    f.eval_fp = &fn1_ginv_eval_fp;
    f.eval_ia = &fn1_ginv_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_ginv_eval_aa;
    f.xd = fn1_ginv_xd;
    f.yd = fn1_ginv_yd;
    f.epsilon = fn1_ginv_epsilon;
    f.delta = fn1_ginv_delta;
    f.nsub = fn1_ginv_nsub;
    return f;
  }
 
