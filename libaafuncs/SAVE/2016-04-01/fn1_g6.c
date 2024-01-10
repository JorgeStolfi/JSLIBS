/* See {fn1_g6.h}. */
/* Last edited on 2005-09-25 15:46:45 by stolfi */
 
#include <fn1_g6.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g6_eval_fp (Float x);
Interval fn1_g6_eval_ia (Interval x);
AAP fn1_g6_eval_aa (AAP x);
 
char *fn1_g6_tag = "g6";
char *fn1_g6_descr = 
  "u = 3(x-1)/4; v = 3x/2;\n"
  "r2 = u^2 + v^2;\n"
  "m2 = (u^2 - 3*v^2)^2;\n"
  "g = u*m2";

Float fn1_g6_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float u = Three*Quarter * (x - One);
      Float v = Three*Half * x;
      Float u2 = u*u;
      Float v2 = v*v;
      Float m2 = u2 - Three * v2;
      Float s = u*m2;
      Float res = s;
      return (res);
    }
  }

Interval fn1_g6_eval_ia (Interval x)
  {
    Interval u = ia_affine(x, Three, Four, -Three*Quarter, Zero);
    Interval v = ia_scale(x, Three, Two);
    Interval u2 = ia_sqr(u);
    Interval v2 = ia_sqr(v);
    Interval m2 = ia_sub(u2, ia_scale(v2, Three, One));
    Interval s = ia_mul(u, m2);
    Interval res = s;
    return (res);
  }

AAP fn1_g6_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP u = aa_affine(x, Three, Four, -Three*Quarter, Zero);
    AAP v = aa_scale(x, Three, Two);
    AAP u2 = aa_sqr(u);
    AAP v2 = aa_sqr(v);
    AAP m2 = aa_sub(u2, aa_scale(v2, Three, One));
    AAP s = aa_mul(u, m2);
    AAP res = s;
    return (aa_return(frame, res));
  }

Interval fn1_g6_xd = {-One, One};
Interval fn1_g6_yd = {-One, One};

Float fn1_g6_epsilon = 1.0e-6;
Float fn1_g6_delta = 1.0e-20;
int fn1_g6_nsub = 16;
 
fn1_data_t fn1_g6_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g6_tag;
    f.descr = fn1_g6_descr;
    f.eval_fp = &fn1_g6_eval_fp;
    f.eval_ia = &fn1_g6_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g6_eval_aa;
    f.xd = fn1_g6_xd;
    f.yd = fn1_g6_yd;
    f.epsilon = fn1_g6_epsilon;
    f.delta = fn1_g6_delta;
    f.nsub = fn1_g6_nsub;
    return f;
  }
 
