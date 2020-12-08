/* See {fn2_f0.h}. */
/* Last edited on 2005-09-25 14:21:57 by stolfi */
 
#include <fn2_f0.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_f0_eval_fp (Float x, Float y);
Interval fn2_f0_eval_ia (Interval x, Interval y);
AAP fn2_f0_eval_aa (AAP x, AAP y);
 
char *fn2_f0_tag = "f0";
char *fn2_f0_descr = "f = x^2 + y^2 + x*y - 1/4";

Float fn2_f0_eval_fp (Float x, Float y)
  { ROUND_NEAR;
    return (x*x + y*y + x*y - Quarter);
  }

Interval fn2_f0_eval_ia (Interval x, Interval y)
  { Interval x2 = ia_sqr(x);
    Interval y2 = ia_sqr(y);
    Interval xy = ia_mul(x, y);
    Interval sum = ia_add(ia_add(x2, y2), xy);
    Interval res = ia_shift (sum, -Quarter);
    return (res);
  }

AAP fn2_f0_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP x2  = aa_sqr(x);
    AAP y2  = aa_sqr(y);
    AAP xy  = aa_mul(x, y);
    AAP sum = aa_add(aa_add(x2, y2), xy);
    AAP res = aa_shift(sum, -Quarter);
    return (aa_return(frame, res));
  }

Interval fn2_f0_xd = {-One, One};
Interval fn2_f0_yd = {-One, One};

int fn2_f0_fn = 32;
 
fn2_data_t fn2_f0_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_f0_tag;
    f.descr = fn2_f0_descr;
    f.eval_fp = &fn2_f0_eval_fp;
    f.eval_ia = &fn2_f0_eval_ia;
    f.eval_aa = &fn2_f0_eval_aa;
    f.xd = fn2_f0_xd;
    f.yd = fn2_f0_yd;
    return f;
  }
 
