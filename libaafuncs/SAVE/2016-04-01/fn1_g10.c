/* See {fn1_g10.h}. */
/* Last edited on 2005-09-25 15:34:15 by stolfi */
 
#include <fn1_g10.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g10_eval_fp (Float x);
Interval fn1_g10_eval_ia (Interval x);
AAP fn1_g10_eval_aa (AAP x);
 
char *fn1_g10_tag = "g10";
char *fn1_g10_descr = "g(x) = 1/sqrt(x^2 + 1/2)";

Float fn1_g10_eval_fp (Float x)
  {
    ROUND_NEAR;
    return (One/sqrt(x*x + Half));
  }

Interval fn1_g10_eval_ia (Interval x)
  {
    Interval x2  = ia_sqr(x);
    Interval hlf = (Interval){Half, Half};
    Interval sum = ia_add(x2, hlf);
    Interval res = ia_inv(ia_sqrt(sum));
    return (res);
  }

AAP fn1_g10_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP x2  = aa_sqr(x);
    AAP hlf = aa_const(Half, Zero);
    AAP sum = aa_add(x2, hlf);
    AAP res = aa_inv(aa_sqrt(sum));
    return (aa_return(frame, res));
  }

Interval fn1_g10_xd = {-Three, Three};
Interval fn1_g10_yd = {-Three, Three};

Float fn1_g10_epsilon = 1.0e-6;
Float fn1_g10_delta = 1.0e-20;
int fn1_g10_nsub = 32;
 
fn1_data_t fn1_g10_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g10_tag;
    f.descr = fn1_g10_descr;
    f.eval_fp = &fn1_g10_eval_fp;
    f.eval_ia = &fn1_g10_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g10_eval_aa;
    f.xd = fn1_g10_xd;
    f.yd = fn1_g10_yd;
    f.epsilon = fn1_g10_epsilon;
    f.delta = fn1_g10_delta;
    f.nsub = fn1_g10_nsub;
    return f;
  }
 
