/* See {fn1_g8.h}. */
/* Last edited on 2005-09-25 15:46:59 by stolfi */
 
#include <fn1_g8.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g8_eval_fp (Float x);
Interval fn1_g8_eval_ia (Interval x);
AAP fn1_g8_eval_aa (AAP x);
 
char *fn1_g8_tag = "g8";
char *fn1_g8_descr = 
  "y = -range(x);\n"
  "d = x-y;\n"
  "t = d(1-d^2);\n"
  "i = 1/(1 + t^2);\n"
  "h = t*i;\n"
  "g = i";

Float fn1_g8_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float y = -x;
      Float d = x - y;
      Float t = Quarter * d * (One - d*d);
      Float i = One / (One + t*t);
      /* Float h = t * i; */
      Float res = i;
      return (res);
    }
  }

Interval fn1_g8_eval_ia (Interval x)
  {
    Interval y = ia_neg(x);
    Interval d = ia_sub(x, y);
    Interval t = ia_mul(d, ia_affine(ia_sqr(d), -One, Four, Quarter, Zero));
    Interval i = ia_inv(ia_shift(ia_sqr(t), One));
    /* Interval h = ia_mul(t, i); */
    Interval res = i;
    return (res);
  }

AAP fn1_g8_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP y = aa_from_interval(ia_neg(aa_range(x)));
    AAP d = aa_sub(x, y);
    AAP t = aa_mul(d, aa_affine(aa_sqr(d), -One, Four, Quarter, Zero));
    AAP i = aa_inv(aa_shift(aa_sqr(t), One));
    /* AAP h = aa_mul(t, i); */
    AAP res = i;
    return (aa_return(frame, res));
  }

Interval fn1_g8_xd = {-Two, Two};
Interval fn1_g8_yd = {-Two, Two};

Float fn1_g8_epsilon = 1.0e-6;
Float fn1_g8_delta = 1.0e-20;
int fn1_g8_nsub = 16;
 
fn1_data_t fn1_g8_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g8_tag;
    f.descr = fn1_g8_descr;
    f.eval_fp = &fn1_g8_eval_fp;
    f.eval_ia = &fn1_g8_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g8_eval_aa;
    f.xd = fn1_g8_xd;
    f.yd = fn1_g8_yd;
    f.epsilon = fn1_g8_epsilon;
    f.delta = fn1_g8_delta;
    f.nsub = fn1_g8_nsub;
    return f;
  }
 
