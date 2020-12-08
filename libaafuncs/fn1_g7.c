/* See {fn1_g7.h}. */
/* Last edited on 2016-12-26 21:34:11 by stolfilocal */
 
#include <fn1_g7.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g7_eval_fp (Float x);
Interval fn1_g7_eval_ia (Interval x);
AAP fn1_g7_eval_aa (AAP x);
 
char *fn1_g7_tag = "g7";
char *fn1_g7_descr = 
  "y = -range(x);\n"
  "d = x-y;\n"
  "t = d(1-d^2);\n"
  "h = t/(1+t^2);\n"
  "u = x+h;\n"
  "v = x+h;\n"
  "g = u^2 + v^2 + 2uv - 1/4";

Float fn1_g7_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float y = -x;
      Float d = x - y;
      Float t = Quarter * d * (One - d*d);
      Float h = t / (1 + t*t);
      Float u = x + h;
      Float v = y + h;
      Float u2 = u*u;
      Float v2 = v*v;
      Float uv = u*v;
      Float m = u2 + v2 + uv + uv;
      Float res = m - Quarter;
      return (res);
    }
  }

Interval fn1_g7_eval_ia (Interval x)
  {
    Interval y = ia_neg(x);
    Interval d = ia_sub(x, y);
    Interval t = ia_mul(d, ia_affine(ia_sqr(d), -One, Four, Quarter, Zero));
    Interval h = ia_div(t, ia_shift(ia_sqr(t), One));
    Interval u = ia_add(x, h);
    Interval v = ia_add(y, h);
    Interval u2 = ia_sqr(u);
    Interval v2 = ia_sqr(v);
    Interval uv = ia_mul(u, v);
    Interval m = ia_add(ia_add(u2, v2), ia_add(uv, uv));
    Interval res = ia_shift(m, - Quarter);
    return (res);
  }

AAP fn1_g7_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP y = aa_from_interval(ia_neg(aa_range(x)));
    AAP d = aa_sub(x, y);
    AAP t = aa_mul(d, aa_affine(aa_sqr(d), -One, Four, Quarter, Zero));
    AAP h = aa_div(t, aa_shift(aa_sqr(t), One));
    AAP u = aa_add(x, h);
    AAP v = aa_add(y, h);
    AAP u2 = aa_sqr(u);
    AAP v2 = aa_sqr(v);
    AAP uv = aa_mul(u, v);
    AAP m = aa_add(aa_add(u2, v2), aa_add(uv, uv));
    AAP res = aa_shift(m, - Quarter);
    return (aa_return(frame, res));
  }

Interval fn1_g7_xd = {-Two, Two};
Interval fn1_g7_yd = {-Two, Two};

Float fn1_g7_epsilon = 1.0e-6f;
Float fn1_g7_delta = 1.0e-20f;
int fn1_g7_nsub = 16;
 
fn1_data_t fn1_g7_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g7_tag;
    f.descr = fn1_g7_descr;
    f.eval_fp = &fn1_g7_eval_fp;
    f.eval_ia = &fn1_g7_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g7_eval_aa;
    f.xd = fn1_g7_xd;
    f.yd = fn1_g7_yd;
    f.epsilon = fn1_g7_epsilon;
    f.delta = fn1_g7_delta;
    f.nsub = fn1_g7_nsub;
    return f;
  }
 
