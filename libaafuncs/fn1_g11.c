/* See {fn1_g11.h}. */
/* Last edited on 2016-12-26 21:34:10 by stolfilocal */
 
#include <fn1_g11.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g11_eval_fp (Float x);
Interval fn1_g11_eval_ia (Interval x);
AAP fn1_g11_eval_aa (AAP x);
 
char *fn1_g11_tag = "g11";
char *fn1_g11_descr = "g(x) = (1 + x)*(1 - x)";

Float fn1_g11_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float u = One + x;
      Float v = One - x;
      Float res = u*v;
      return (res);
    }
  }

Interval fn1_g11_eval_ia (Interval x)
  {
    Interval u = ia_shift(x, One);
    Interval v = ia_shift(ia_neg(x), One);
    Interval res = ia_mul(u, v);
    return (res);
  }

AAP fn1_g11_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP u = aa_shift(x, One);
    AAP v = aa_shift(aa_neg(x), One);
    AAP res = aa_mul(u, v);
    return (aa_return(frame, res));
  }

Interval fn1_g11_xd = {-Two, Two};
Interval fn1_g11_yd = {-Two, Two};

Float fn1_g11_epsilon = 1.0e-6f;
Float fn1_g11_delta = 1.0e-20f;
int fn1_g11_nsub = 7;
 
fn1_data_t fn1_g11_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g11_tag;
    f.descr = fn1_g11_descr;
    f.eval_fp = &fn1_g11_eval_fp;
    f.eval_ia = &fn1_g11_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g11_eval_aa;
    f.xd = fn1_g11_xd;
    f.yd = fn1_g11_yd;
    f.epsilon = fn1_g11_epsilon;
    f.delta = fn1_g11_delta;
    f.nsub = fn1_g11_nsub;
    return f;
  }
 
