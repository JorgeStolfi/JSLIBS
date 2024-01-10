/* See {fn2_f5.h}. */
/* Last edited on 2016-12-26 21:27:35 by stolfilocal */
 
#include <fn2_f5.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_f5_eval_fp (Float x, Float y);
Interval fn2_f5_eval_ia (Interval x, Interval y);
AAP fn2_f5_eval_aa (AAP x, AAP y);
 
char *fn2_f5_tag = "f5";
char *fn2_f5_descr = 
  "r2 = x^2 + y^2;\nm2 = (x^2 - 3*y^2)^2;\ns = x*m2 - r2^2;\nt = (1/8)/(2 + r2);\nf = s + t";
 
Float fn2_f5_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float x2 = x*x;
      Float y2 = y*y;
      Float r2 = x2 + y2;
      Float m2 = x2 - Three * y2;
      Float s = x*m2 + r2*r2;
      Float t = (Float)(1.0/(8.0*(2.0 + r2)));
      Float res = s + t;
      return (res);
    }
  }

Interval fn2_f5_eval_ia (Interval x, Interval y)
  {
    Interval x2 = ia_sqr(x);
    Interval y2 = ia_sqr(y);
    Interval r2 = ia_add(x2, y2);
    Interval m2 = ia_sub(x2, ia_scale(y2, Three, One));
    Interval s = ia_add(ia_mul(x, m2), ia_sqr(r2));
    Interval t = ia_scale(ia_inv(ia_shift(r2, Two)), One, Eight);
    Interval res = ia_add(s, t);
    return (res);
  }

AAP fn2_f5_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP x2 = aa_sqr(x);
    AAP y2 = aa_sqr(y);
    AAP r2 = aa_add(x2, y2);
    AAP m2 = aa_sub(x2, aa_scale(y2, Three, One));
    AAP s = aa_add(aa_mul(x, m2), aa_sqr(r2));
    AAP t = aa_scale(aa_inv(aa_shift(r2, Two)), One, Eight);
    AAP res = aa_add(s, t);
    return (aa_return(frame, res));
  }

Interval fn2_f5_xd = {-Three/Two, Three/Two};
Interval fn2_f5_yd = {-Three/Two, Three/Two};

int fn2_f5_fn = 32;
 
fn2_data_t fn2_f5_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_f5_tag;
    f.descr = fn2_f5_descr;
    f.eval_fp = &fn2_f5_eval_fp;
    f.eval_ia = &fn2_f5_eval_ia;
    f.eval_aa = &fn2_f5_eval_aa;
    f.xd = fn2_f5_xd;
    f.yd = fn2_f5_yd;
    return f;
  }
 
