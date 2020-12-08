/* See {fn1_g12.h}. */
/* Last edited on 2016-12-26 21:34:10 by stolfilocal */
 
#include <fn1_g12.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g12_eval_fp (Float x);
Interval fn1_g12_eval_ia (Interval x);
AAP fn1_g12_eval_aa (AAP x);
 
char *fn1_g12_tag = "g12";
char *fn1_g12_descr = 
  "y = [x]+2;\n"
  "d = x-y;\n"
  "t = d(1-d^2)/4; h = t/(1+t^2);\n"
  "u = x+h; v = x+h;\n"
  "m = u^2 + v^2 + 2uv;\n"
  "g = m - 1/4";

Float fn1_g12_eval_fp_2 (Float x, Float y);
Interval fn1_g12_eval_ia_2 (Interval x, Interval y);
AAP fn1_g12_eval_aa_2 (AAP x, AAP y);

Float fn1_g12_eval_fp_2 (Float x, Float y)
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
  
Float fn1_g12_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float y = x + Two;
      Float res = fn1_g12_eval_fp_2(x, y);
      return (res);
    }
  }
    
Interval fn1_g12_eval_ia_2 (Interval x, Interval y)
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
    Interval res = ia_shift(m, -Quarter);
    return (res);
  }
  
Interval fn1_g12_eval_ia (Interval x)
  {
    Interval y = ia_shift(x, Two);
    Interval res = fn1_g12_eval_ia_2(x, y);
    return (res);
  }

AAP fn1_g12_eval_aa_2 (AAP x, AAP y)
  {
    AAP d = aa_sub(x, y);
    AAP t = aa_mul(d, aa_affine(aa_sqr(d), -One, Four, Quarter, Zero));
    AAP h = aa_div(t, aa_shift(aa_sqr(t), One));
    AAP u = aa_add(x, h);
    AAP v = aa_add(y, h);
    AAP u2 = aa_sqr(u);
    AAP v2 = aa_sqr(v);
    AAP uv = aa_mul(u, v);
    AAP m = aa_add(aa_add(u2, v2), aa_add(uv, uv));
    AAP res = aa_shift(m, -Quarter);
    return (res);
  }

AAP fn1_g12_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP y = aa_shift(aa_from_interval(aa_range(x)), Two);
    AAP res = fn1_g12_eval_aa_2(x, y);
    return (aa_return(frame, res));
  }

Interval fn1_g12_xd = {-Four, Four};
Interval fn1_g12_yd = {-3.0, 15.0};

Float fn1_g12_epsilon = 1.0e-6f;
Float fn1_g12_delta = 1.0e-20f;
int fn1_g12_nsub = 8;
 
fn1_data_t fn1_g12_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g12_tag;
    f.descr = fn1_g12_descr;
    f.eval_fp = &fn1_g12_eval_fp;
    f.eval_ia = &fn1_g12_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g12_eval_aa;
    f.xd = fn1_g12_xd;
    f.yd = fn1_g12_yd;
    f.epsilon = fn1_g12_epsilon;
    f.delta = fn1_g12_delta;
    f.nsub = fn1_g12_nsub;
    return f;
  }
 
