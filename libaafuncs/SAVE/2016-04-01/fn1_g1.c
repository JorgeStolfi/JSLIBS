/* See {fn1_g1.h}. */
/* Last edited on 2005-09-25 15:52:44 by stolfi */
 
#include <fn1_g1.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g1_eval_fp (Float x);
Interval fn1_g1_eval_ia (Interval x);
AAP fn1_g1_eval_aa (AAP x);
 
char *fn1_g1_tag = "g1";
char *fn1_g1_descr = "g(x) = sqrt(x^2 - x + 1/2)/sqrt(x^2 + 1/2)";

Float fn1_g1_eval_fp (Float x)
  {
    ROUND_NEAR;
    return (sqrt(x*x - x + 0.50)/sqrt(x*x + 0.50));
  }

Interval fn1_g1_eval_ia (Interval x)
  {
    Interval x2  = ia_sqr(x);
    Interval hlf = {Half, Half};
    Interval dif = ia_sub(x2, x);
    Interval sum1 = ia_add(dif, hlf);
    Interval sum2 = ia_add(x2, hlf);
    Interval res = ia_div(ia_sqrt(sum1), ia_sqrt(sum2));
    return (res);
  }

AAP fn1_g1_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP x2  = aa_sqr(x);
    AAP hlf = aa_const(Half, Zero);
    AAP dif = aa_sub(x2, x);
    AAP sum1 = aa_add(dif, hlf);
    AAP sum2 = aa_add(x2, hlf);
    AAP res = aa_div(aa_sqrt(sum1), aa_sqrt(sum2));
    return (aa_return(frame, res));
  }

Interval fn1_g1_xd = {-Two, Two};
Interval fn1_g1_yd = {-Quarter, Two};

Float fn1_g1_epsilon = 1.0e-6;
Float fn1_g1_delta = 1.0e-20;
int fn1_g1_nsub = 16;
 
fn1_data_t fn1_g1_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g1_tag;
    f.descr = fn1_g1_descr;
    f.eval_fp = &fn1_g1_eval_fp;
    f.eval_ia = &fn1_g1_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g1_eval_aa;
    f.xd = fn1_g1_xd;
    f.yd = fn1_g1_yd;
    f.epsilon = fn1_g1_epsilon;
    f.delta = fn1_g1_delta;
    f.nsub = fn1_g1_nsub;
    return f;
  }
 
