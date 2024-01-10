/* See {fn2_f2.h}. */
/* Last edited on 2005-09-25 14:22:11 by stolfi */
 
#include <fn2_f2.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_f2_eval_fp (Float x, Float y);
Interval fn2_f2_eval_ia (Interval x, Interval y);
AAP fn2_f2_eval_aa (AAP x, AAP y);
 
char *fn2_f2_tag = "f2";
char *fn2_f2_descr = 
  "r2 = x^2 + y^2;\nxy = x*y;\nfa = (r^2 + xy - 1/4)^3;\nfb = (r^2 - xy - 1/4)^3;\nf = fa + fb - fa*fb";

Float fn2_f2_eval_fp (Float x, Float y)
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

Interval fn2_f2_eval_ia (Interval x, Interval y)
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

AAP fn2_f2_eval_aa (AAP x, AAP y)
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

Interval fn2_f2_xd = {-One, One};
Interval fn2_f2_yd = {-One, One};

int fn2_f2_fn = 32;
 
fn2_data_t fn2_f2_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_f2_tag;
    f.descr = fn2_f2_descr;
    f.eval_fp = &fn2_f2_eval_fp;
    f.eval_ia = &fn2_f2_eval_ia;
    f.eval_aa = &fn2_f2_eval_aa;
    f.xd = fn2_f2_xd;
    f.yd = fn2_f2_yd;
    return f;
  }
 
