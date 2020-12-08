/* See {fn1_gbadia.h}. */
/* Last edited on 2016-12-26 21:27:51 by stolfilocal */
 
#include <fn1_gbadia.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_gbadia_eval_fp (Float x);
Interval fn1_gbadia_eval_ia (Interval x);
AAP fn1_gbadia_eval_aa (AAP x);
 
char *fn1_gbadia_tag = "gbadia";
char *fn1_gbadia_descr = "g(x) = sqrt(x^2 - x + 1/2)/sqrt(x^2 + 1/2)";

Float fn1_gbadia_eval_fp (Float x)
  {
    ROUND_NEAR;
    return (Float)(sqrt(x*x - x + 0.50)/sqrt(x*x + 0.50));
  }

Interval fn1_gbadia_eval_ia (Interval x)
  {
    Interval x2  = ia_sqr(x);
    Interval hlf = {Half, Half};
    Interval dif = ia_sub(x2, x);
    Interval sum1 = ia_add(dif, hlf);
    Interval sum2 = ia_add(x2, hlf);
    Interval res = ia_div(ia_sqrt(sum1), ia_sqrt(sum2));
    return (res);
  }

AAP fn1_gbadia_eval_aa (AAP x)
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

Interval fn1_gbadia_xd = {-Two, Two};
Interval fn1_gbadia_yd = {-Half*Quarter, Two};

Float fn1_gbadia_epsilon = 1.0e-6f;
Float fn1_gbadia_delta = 1.0e-20f;
int fn1_gbadia_nsub = 16;
 
fn1_data_t fn1_gbadia_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_gbadia_tag;
    f.descr = fn1_gbadia_descr;
    f.eval_fp = &fn1_gbadia_eval_fp;
    f.eval_ia = &fn1_gbadia_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_gbadia_eval_aa;
    f.xd = fn1_gbadia_xd;
    f.yd = fn1_gbadia_yd;
    f.epsilon = fn1_gbadia_epsilon;
    f.delta = fn1_gbadia_delta;
    f.nsub = fn1_gbadia_nsub;
    return f;
  }
 
