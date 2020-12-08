/* See {fn2_f4.h}. */
/* Last edited on 2005-09-25 14:22:22 by stolfi */
 
#include <fn2_f4.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_f4_eval_fp (Float x, Float y);
Interval fn2_f4_eval_ia (Interval x, Interval y);
AAP fn2_f4_eval_aa (AAP x, AAP y);
 
char *fn2_f4_tag = "f4";
char *fn2_f4_descr = 
  "r2 = x^2 + y^2;\nm2 = (x^2 - 3*y^2)^2;\nf = x*m2 - r2^2";

Float fn2_f4_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float x2 = x*x;
      Float y2 = y*y;
      Float r2 = x2 + y2;
      Float m2 = x2 - Three * y2;
      Float res = x*m2 + r2*r2;
      return (res);
    }
  }

Interval fn2_f4_eval_ia (Interval x, Interval y)
  {
    Interval x2 = ia_sqr(x);
    Interval y2 = ia_sqr(y);
    Interval r2 = ia_add(x2, y2);
    Interval m2 = ia_sub(x2, ia_scale(y2, Three, One));
    Interval res = ia_add(ia_mul(x, m2), ia_sqr(r2));
    return (res);
  }

AAP fn2_f4_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP x2 = aa_sqr(x);
    AAP y2 = aa_sqr(y);
    AAP r2 = aa_add(x2, y2);
    AAP m2 = aa_sub(x2, aa_scale(y2, Three, One));
    AAP res = aa_add(aa_mul(x, m2), aa_sqr(r2));
    return (aa_return(frame, res));
  }

Interval fn2_f4_xd = {-Three/Two, Three/Two};
Interval fn2_f4_yd = {-Three/Two, Three/Two};

int fn2_f4_fn = 32;
 
fn2_data_t fn2_f4_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_f4_tag;
    f.descr = fn2_f4_descr;
    f.eval_fp = &fn2_f4_eval_fp;
    f.eval_ia = &fn2_f4_eval_ia;
    f.eval_aa = &fn2_f4_eval_aa;
    f.xd = fn2_f4_xd;
    f.yd = fn2_f4_yd;
    return f;
  }
 
