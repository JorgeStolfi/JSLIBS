/* See {fn2_f7.h}. */
/* Last edited on 2005-09-25 14:22:37 by stolfi */
 
#include <fn2_f7.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_f7_eval_fp (Float x, Float y);
Interval fn2_f7_eval_ia (Interval x, Interval y);
AAP fn2_f7_eval_aa (AAP x, AAP y);
 
char *fn2_f7_tag = "f7";
char *fn2_f7_descr = 
  "d = x-y; t = d(1-d^2); h = t/(1+t^2); u = x+h; v = x+h; f = u^2 + v^2 + 2uv - 1/4";

Float fn2_f7_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    {
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

Interval fn2_f7_eval_ia (Interval x, Interval y)
  {
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

AAP fn2_f7_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
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

Interval fn2_f7_xd = {-Four, Four};
Interval fn2_f7_yd = {-Four, Four};

int fn2_f7_fn = 32;
 
fn2_data_t fn2_f7_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_f7_tag;
    f.descr = fn2_f7_descr;
    f.eval_fp = &fn2_f7_eval_fp;
    f.eval_ia = &fn2_f7_eval_ia;
    f.eval_aa = &fn2_f7_eval_aa;
    f.xd = fn2_f7_xd;
    f.yd = fn2_f7_yd;
    return f;
  }
 
