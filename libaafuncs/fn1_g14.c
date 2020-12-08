/* See {fn1_g14.h}. */
/* Last edited on 2016-12-26 21:34:07 by stolfilocal */
 
#include <fn1_g14.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g14_eval_fp (Float x);
Interval fn1_g14_eval_ia (Interval x);
AAP fn1_g14_eval_aa (AAP x);
 
char *fn1_g14_tag = "g14";
char *fn1_g14_descr = 
  "y = 1/32 +/- 1/32;\n"
  "r2 = x^2 + y^2;\n"
  "xy = x*y;\n"
  "fa = (r^2 + xy - 1/4)^3;\n"
  "fb = (r^2 - xy - 1/4)^3;\n"
  "g = fa + fb - fa*fb";

Float fn1_g14_eval_fp_2 (Float x, Float y);
Interval fn1_g14_eval_ia_2 (Interval x, Interval y);
AAP fn1_g14_eval_aa_2 (AAP x, AAP y);

Float fn1_g14_eval_fp_2 (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float ha = x*x + y*y + x*y - Quarter;
      Float fa = ha*ha*ha;
      Float hb = x*x + y*y - x*y - Quarter;
      Float fb = hb*hb*hb;
      Float res = (fa + fb - fa*fb);
      return (res);
    }
  }

Float fn1_g14_eval_fp (Float x)
  {
    Float y = (Float)((1.0 + sin(32.0*M_PI*x))/32.0);
    return(fn1_g14_eval_fp_2(x, y));
  }

Interval fn1_g14_eval_ia_2 (Interval x, Interval y)
  {
    Interval x2 = ia_sqr(x);
    Interval y2 = ia_sqr(y);
    Interval r2 = ia_add(x2, y2);
    Interval xy = ia_mul(x, y);
    Interval ha = ia_shift(ia_add(r2, xy), -Quarter);
    Interval fa = ia_mul(ia_sqr(ha), ha);
    Interval hb = ia_shift(ia_sub(r2, xy), -Quarter);
    Interval fb = ia_mul(ia_sqr(hb), hb);
    Interval fafb = ia_mul(fa, fb);
    Interval res = ia_sub (ia_add(fa, fb), fafb);
    return (res);
  }
  
Interval fn1_g14_eval_ia (Interval x)
  {
    Interval y = ia_const(1/32.0, 1/32.0);
    return(fn1_g14_eval_ia_2(x, y));
  }

AAP fn1_g14_eval_aa_2 (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP x2 = aa_sqr(x);
    AAP y2 = aa_sqr(y);
    AAP r2 = aa_add(x2, y2);
    AAP xy = aa_mul(x, y);
    AAP ha = aa_shift(aa_add(r2, xy), -Quarter);
    AAP fa = aa_mul(aa_sqr(ha), ha);
    AAP hb = aa_shift(aa_sub(r2, xy), -Quarter);
    AAP fb = aa_mul(aa_sqr(hb), hb);
    AAP fafb = aa_mul(fa, fb);
    AAP res = aa_sub (aa_add(fa, fb), fafb);
    return (aa_return(frame, res));
  }

AAP fn1_g14_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP y = aa_const(1/32.0, 1/32.0);
    return(aa_return(frame, fn1_g14_eval_aa_2(x, y)));
  }

Interval fn1_g14_xd = {-One, One};
Interval fn1_g14_yd = {-One, One};

Float fn1_g14_epsilon = 1.0e-6f;
Float fn1_g14_delta = 1.0e-20f;
int fn1_g14_nsub = 32;
 
fn1_data_t fn1_g14_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g14_tag;
    f.descr = fn1_g14_descr;
    f.eval_fp = &fn1_g14_eval_fp;
    f.eval_ia = &fn1_g14_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g14_eval_aa;
    f.xd = fn1_g14_xd;
    f.yd = fn1_g14_yd;
    f.epsilon = fn1_g14_epsilon;
    f.delta = fn1_g14_delta;
    f.nsub = fn1_g14_nsub;
    return f;
  }
 
