/* See {fn1_g5.h}. */
/* Last edited on 2016-12-26 21:34:32 by stolfilocal */
 
#include <fn1_g5.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g5_eval_fp (Float x);
Interval fn1_g5_eval_ia (Interval x);
AAP fn1_g5_eval_aa (AAP x);
 
char *fn1_g5_tag = "g5";
char *fn1_g5_descr = 
  "u = 3(x-1)/4; v = 3x/2; r2 = u^2 + v^2;\n"
  "m2 = (u^2 - 3*v^2)^2;\n"
  "s = u*m2 - r2^2;\n"
  "t = (1/8)/(2 + r2);\n"
  "g = s + t";

Float fn1_g5_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float u = Three*Quarter * (x - One);
      Float v = Three*Half * x;
      Float u2 = u*u;
      Float v2 = v*v;
      Float r2 = u2 + v2;
      Float m2 = u2 - Three * v2;
      Float s = u*m2 + r2*r2;
      Float t = (Float)((1.0/8.0)/(2.0 + r2));
      Float res = s + t;
      return (res);
    }
  }

Interval fn1_g5_eval_ia (Interval x)
  {
    Float Eighth = Half*Quarter;
    Interval u = ia_affine(x, Three, Four, -Three*Quarter, Zero);
    Interval v = ia_scale(x, Three, Two);
    Interval u2 = ia_sqr(u);
    Interval v2 = ia_sqr(v);
    Interval r2 = ia_add(u2, v2);
    Interval m2 = ia_sub(u2, ia_scale(v2, Three, One));
    Interval s = ia_add(ia_mul(u, m2), ia_sqr(r2));
    Interval t = ia_scale(ia_inv(ia_shift(r2, Two)), Eighth, One);
    Interval res = ia_add(s, t);
    return (res);
  }

AAP fn1_g5_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    Float Eighth = Half*Quarter;
    AAP u = aa_affine(x, Three, Four, -Three*Quarter, Zero);
    AAP v = aa_scale(x, Three, Two);
    AAP u2 = aa_sqr(u);
    AAP v2 = aa_sqr(v);
    AAP r2 = aa_add(u2, v2);
    AAP m2 = aa_sub(u2, aa_scale(v2, Three, One));
    AAP s = aa_add(aa_mul(u, m2), aa_sqr(r2));
    AAP t = aa_scale(aa_inv(aa_shift(r2, Two)), Eighth, One);
    AAP res = aa_add(s, t);
    return (aa_return(frame, res));
  }

Interval fn1_g5_xd = {-One, One};
Interval fn1_g5_yd = {-One, One};

Float fn1_g5_epsilon = 1.0e-6f;
Float fn1_g5_delta = 1.0e-20f;
int fn1_g5_nsub = 16;
 
fn1_data_t fn1_g5_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g5_tag;
    f.descr = fn1_g5_descr;
    f.eval_fp = &fn1_g5_eval_fp;
    f.eval_ia = &fn1_g5_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g5_eval_aa;
    f.xd = fn1_g5_xd;
    f.yd = fn1_g5_yd;
    f.epsilon = fn1_g5_epsilon;
    f.delta = fn1_g5_delta;
    f.nsub = fn1_g5_nsub;
    return f;
  }
 
