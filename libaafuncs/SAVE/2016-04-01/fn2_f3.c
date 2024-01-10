/* See {fn2_f3.h}. */
/* Last edited on 2005-09-25 14:22:16 by stolfi */
 
#include <fn2_f3.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_f3_eval_fp (Float x, Float y);
Interval fn2_f3_eval_ia (Interval x, Interval y);
AAP fn2_f3_eval_aa (AAP x, AAP y);
 
char *fn2_f3_tag = "f3";
char *fn2_f3_descr = 
  "r2 = x^2 + y^2;\nm2 = (x^2 - 3*y^2)^2;\nf = x^2*m^2 - r^2 + 1/4";

Float fn2_f3_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float x2 = x*x;
      Float y2 = y*y;
      Float r2 = x2 + y2;
      Float m = x2 - Three * y2;
      Float m2 = m*m;
      Float res = x2*m2 - r2 + Quarter;
      return (res);
    }
  }

Interval fn2_f3_eval_ia (Interval x, Interval y)
  {
    Interval x2 = ia_sqr(x);
    Interval y2 = ia_sqr(y);
    Interval r2 = ia_add(x2, y2);
    Interval m = ia_sub(x2, ia_scale(y2, Three, One));
    Interval m2 = ia_sqr(m);
    Interval res = ia_shift(ia_sub(ia_mul(x2, m2),  r2), Quarter);
    return (res);
  }

AAP fn2_f3_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP x2 = aa_sqr(x);
    AAP y2 = aa_sqr(y);
    AAP r2 = aa_add(x2, y2);
    AAP m = aa_sub(x2, aa_scale(y2, Three, One));
    AAP m2 = aa_sqr(m);
    AAP res = aa_shift(aa_sub(aa_mul(x2, m2),  r2), Quarter);
    return (aa_return(frame, res));
  }

Interval fn2_f3_xd = {-Three, Three};
Interval fn2_f3_yd = {-Three, Three};

int fn2_f3_fn = 64;
 
fn2_data_t fn2_f3_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_f3_tag;
    f.descr = fn2_f3_descr;
    f.eval_fp = &fn2_f3_eval_fp;
    f.eval_ia = &fn2_f3_eval_ia;
    f.eval_aa = &fn2_f3_eval_aa;
    f.xd = fn2_f3_xd;
    f.yd = fn2_f3_yd;
    return f;
  }
 
