/* See {fn1_gbadmul.h}. */
/* Last edited on 2016-12-26 21:33:28 by stolfilocal */
 
#include <fn1_gbadmul.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_gbadmul_eval_fp (Float x);
Interval fn1_gbadmul_eval_ia (Interval x);
AAP fn1_gbadmul_eval_aa (AAP x);
 
char *fn1_gbadmul_tag = "gbadmul";
char *fn1_gbadmul_descr = 
  "g(x) = (10 + x)*(10 - x)";

Float fn1_gbadmul_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float u = 10.0f + x;
      Float v = 10.0f - x;
      Float res = u*v/10.0f;
      return (res);
    }
  }

Interval fn1_gbadmul_eval_ia (Interval x)
  {
    Interval u = ia_shift(x, 10.0);
    Interval v = ia_shift(ia_neg(x), 10.0);
    Interval res = ia_scale(ia_mul(u, v), One, 10.0);
    return (res);
  }

AAP fn1_gbadmul_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP u = aa_shift(x, 10.0);
    AAP v = aa_shift(aa_neg(x), 10.0);
    AAP res = aa_scale(aa_mul(u, v), One, 10.0);
    return (aa_return(frame, res));
  }

Interval fn1_gbadmul_xd = {-22.0, 22.0};
Interval fn1_gbadmul_yd = {-12.0, 16.0};

Float fn1_gbadmul_epsilon = 1.0e-6f;
Float fn1_gbadmul_delta = 1.0e-20f;
int fn1_gbadmul_nsub = 11;
 
fn1_data_t fn1_gbadmul_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_gbadmul_tag;
    f.descr = fn1_gbadmul_descr;
    f.eval_fp = &fn1_gbadmul_eval_fp;
    f.eval_ia = &fn1_gbadmul_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_gbadmul_eval_aa;
    f.xd = fn1_gbadmul_xd;
    f.yd = fn1_gbadmul_yd;
    f.epsilon = fn1_gbadmul_epsilon;
    f.delta = fn1_gbadmul_delta;
    f.nsub = fn1_gbadmul_nsub;
    return f;
  }
 
