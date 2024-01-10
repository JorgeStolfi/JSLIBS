/* See {fn1_giaboom.h}. */
/* Last edited on 2005-09-25 15:47:38 by stolfi */
 
#include <fn1_giaboom.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_giaboom_eval_fp (Float x);
Interval fn1_giaboom_eval_ia (Interval x);
AAP fn1_giaboom_eval_aa (AAP x);
 
char *fn1_giaboom_tag = "giaboom";
char *fn1_giaboom_descr = "h(x) = sqrt(x^2 - x + 1/2)/sqrt(x^2 + 1/2) g(x) = h(h(x))";

Float aat_h_fp (Float x);
Interval aat_h_ia (Interval x);
AAP aat_h_aa (AAP x);

Float aat_h_fp (Float x)
  {
    ROUND_NEAR;
    return (sqrt(x*x - x + 0.50)/sqrt(x*x + 0.50));
  }

Float fn1_giaboom_eval_fp (Float x)
  {
    ROUND_NEAR;
    return aat_h_fp(aat_h_fp(x));
  }

Interval aat_h_ia (Interval x)
  {
    Interval x2  = ia_sqr(x);
    Interval hlf = {Half, Half};
    Interval dif = ia_sub(x2, x);
    Interval sum1 = ia_add(dif, hlf);
    Interval sum2 = ia_add(x2, hlf);
    Interval res = ia_div(ia_sqrt(sum1), ia_sqrt(sum2));
    return (res);
  }

Interval fn1_giaboom_eval_ia (Interval x)
  {
    return (aat_h_ia(aat_h_ia(x)));
  }

AAP aat_h_aa (AAP x)
  {
    AAP x2  = aa_sqr(x);
    AAP hlf = aa_const(Half, Zero);
    AAP dif = aa_sub(x2, x);
    AAP sum1 = aa_add(dif, hlf);
    AAP sum2 = aa_add(x2, hlf);
    AAP res = aa_div(aa_sqrt(sum1), aa_sqrt(sum2));
    return (res);
  }

AAP fn1_giaboom_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP res = aat_h_aa(aat_h_aa(x));
    return (aa_return(frame, res));
  }

Interval fn1_giaboom_xd = {-Two, Two};
Interval fn1_giaboom_yd = {-Half*Quarter, Two};

Float fn1_giaboom_epsilon = 1.0e-6;
Float fn1_giaboom_delta = 1.0e-20;
int fn1_giaboom_nsub = 16;
 
fn1_data_t fn1_giaboom_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_giaboom_tag;
    f.descr = fn1_giaboom_descr;
    f.eval_fp = &fn1_giaboom_eval_fp;
    f.eval_ia = &fn1_giaboom_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_giaboom_eval_aa;
    f.xd = fn1_giaboom_xd;
    f.yd = fn1_giaboom_yd;
    f.epsilon = fn1_giaboom_epsilon;
    f.delta = fn1_giaboom_delta;
    f.nsub = fn1_giaboom_nsub;
    return f;
  }
 
