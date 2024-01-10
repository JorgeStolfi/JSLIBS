/* See {fn1_gsin3.h}. */
/* Last edited on 2005-09-25 21:19:40 by stolfi */
 
#include <fn1_gsin3.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_gsin3_eval_fp (Float x);
Interval fn1_gsin3_eval_ia (Interval x);
Interval fn1_gsin3_diff_ia (Interval x);
AAP fn1_gsin3_eval_aa (AAP x);
 
char *fn1_gsin3_tag = "gsin3";
char *fn1_gsin3_descr = 
  "x2 = x^2;\n"
  "t1 = x;\n"
  "t3 = -t1*x2/6;\n"
  "t5 = -t3*x2/20;\n"
  "g(x) = t1 + t2 + t3";

Float fn1_gsin3_eval_fp (Float x)
  {
    ROUND_NEAR;
    { 
      Float x2 = x*x;
      Float t1 = x;
      Float t3 = -t1*x2/6;
      Float t5 = -t3*x2/20;
      Float f = t1 + t3 + t5;
      return f;
    }
  }

Interval fn1_gsin3_eval_ia (Interval x)
  {
    Interval x2 = ia_sqr(x);
    Interval t1 = x;
    Interval t3 = ia_scale(ia_mul(t1, x2), 1, -6);
    Interval t5 = ia_scale(ia_mul(t3, x2), 1, -20);
    Interval f = ia_add(ia_add(t1, t3), t5);
    return f;
  }

Interval fn1_gsin3_diff_ia (Interval x)
  {
    Interval x2 = ia_sqr(x);
    Interval d_t1 = ia_const(1, 0);
    Interval d_t3 = ia_scale(x2, 1, -2);
    Interval d_t5 = ia_scale(d_t3, 1, -12);
    Interval f = ia_add(ia_add(d_t1, d_t3), d_t5);
    return f;
  }

AAP fn1_gsin3_eval_aa (AAP x)
  { 
    MemP frame = aa_top();
    AAP x2 = aa_sqr(x);
    AAP t1 = x;
    AAP t3 = aa_affine(aa_mul(t1, x2), 1, -6, 0, 0);
    AAP t5 = aa_affine(aa_mul(t3, x2), 1, -20, 0, 0);
    AAP f = aa_add(aa_add(t1, t3), t5);
    return aa_return(frame, f);
  }

Interval fn1_gsin3_xd = {-2, 2};
Interval fn1_gsin3_yd = {-2, 2};

Float fn1_gsin3_epsilon = 1.0e-6;
Float fn1_gsin3_delta = 1.0e-20;
int fn1_gsin3_nsub = 16;
 
fn1_data_t fn1_gsin3_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_gsin3_tag;
    f.descr = fn1_gsin3_descr;
    f.eval_fp = &fn1_gsin3_eval_fp;
    f.eval_ia = &fn1_gsin3_eval_ia;
    f.diff_ia = &fn1_gsin3_diff_ia;
    f.eval_aa = &fn1_gsin3_eval_aa;
    f.xd = fn1_gsin3_xd;
    f.yd = fn1_gsin3_yd;
    f.epsilon = fn1_gsin3_epsilon;
    f.delta = fn1_gsin3_delta;
    f.nsub = fn1_gsin3_nsub;
    return f;
  }
 
