/* See {fn1_g2.h}. */
/* Last edited on 2016-12-26 21:34:08 by stolfilocal */
 
#include <fn1_g2.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g2_eval_fp (Float x);
Interval fn1_g2_eval_ia (Interval x);
AAP fn1_g2_eval_aa (AAP x);
 
char *fn1_g2_tag = "g2";
char *fn1_g2_descr = 
  "y = 27/16;\n"
  "r2 = x^2 + y^2;\n"
  "m2 = (x^2 - 3*y^2)^2;\n"
  "g = (x^2*m2 - r2)/20";

#define Twenty (20.0f)

Float fn1_g2_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float yc = 27.0/16.0; /* Should be exact */
      Float y = yc;
      Float x2 = x*x;
      Float y2 = y*y;
      Float r2 = x2 + y2;
      Float m = x2 - Three * y2;
      Float m2 = m*m;
      Float res = (x2*m2 - r2)/Twenty;
      return (res);
    }
  }

Interval fn1_g2_eval_ia (Interval x)
  {
    ROUND_NEAR;
    {
      Float yc = 27.0/16.0; /* Should be exact */
      Interval y = {yc, yc};
      Interval x2 = ia_sqr(x);
      Interval y2 = ia_sqr(y);
      Interval r2 = ia_add(x2, y2);
      Interval m = ia_sub(x2, ia_scale(y2, Three, One));
      Interval m2 = ia_sqr(m);
      Interval res = ia_scale(ia_sub(ia_mul(x2, m2),  r2), One, Twenty);
      return (res);
    }
  }

AAP fn1_g2_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    ROUND_NEAR;
    {
      Float yc = 27.0/16.0; /* Should be exact */
      AAP y = aa_const(yc, Zero);
      AAP x2 = aa_sqr(x);
      AAP y2 = aa_sqr(y);
      AAP r2 = aa_add(x2, y2);
      AAP m = aa_sub(x2, aa_scale(y2, Three, One));
      AAP m2 = aa_sqr(m);
      AAP res = aa_scale(aa_sub(aa_mul(x2, m2),  r2), One, Twenty);
      return (aa_return(frame, res));
    }
  }

Interval fn1_g2_xd = {-6.0, 6.0};
Interval fn1_g2_yd = {-6.0, 6.0};

Float fn1_g2_epsilon = 1.0e-6f;
Float fn1_g2_delta = 1.0e-20f;
int fn1_g2_nsub = 64;
 
fn1_data_t fn1_g2_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g2_tag;
    f.descr = fn1_g2_descr;
    f.eval_fp = &fn1_g2_eval_fp;
    f.eval_ia = &fn1_g2_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g2_eval_aa;
    f.xd = fn1_g2_xd;
    f.yd = fn1_g2_yd;
    f.epsilon = fn1_g2_epsilon;
    f.delta = fn1_g2_delta;
    f.nsub = fn1_g2_nsub;
    return f;
  }
 
