/* See {fn1_f2.h}. */
/* Last edited on 2016-12-26 21:31:57 by stolfilocal */
 
#include <fn1_f2.h>
#include <fn1_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
 
Float fn1_f2_eval_fp (Float x);
Interval fn1_f2_eval_ia (Interval x);
Interval fn1_f2_diff_ia (Interval x);
AAP fn1_f2_eval_aa (AAP x);
 
char *fn1_f2_tag = "f2";
char *fn1_f2_descr = 
  "f(x) = (x - 2)(x + 2)/4";

Float fn1_f2_eval_fp (Float x)
  {
    ROUND_NEAR;
    return (Float)((x-2.0)*(x+2.0)/4.0);
  }

Interval fn1_f2_eval_ia (Interval x)
  {
    Interval f1 = ia_shift(x, -Two);
    Interval f2 = ia_shift(x, +Two);
    Interval f = ia_scale(ia_mul(f1, f2), One, Four);
    return f;
  }

Interval fn1_f2_diff_ia (Interval x)
  {
    Interval one = {One, One};

    Interval f1 = ia_shift(x, -Two);
    Interval d_f1 = one;

    Interval f2 = ia_shift(x, +Two);
    Interval d_f2 = one;

    Interval d_f = ia_scale(ia_add(ia_mul(f1, d_f2), ia_mul(d_f1, f2)), One, Four);

    return d_f;
  }

AAP fn1_f2_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP f1 = aa_shift(x, -Two);
    AAP f2 = aa_shift(x, +Two);
    AAP f = aa_scale(aa_mul(f1, f2), One, Four);
    AAP res = f;
    return (aa_return(frame, res));
  }

Interval fn1_f2_xd = {-Four, Three};
Interval fn1_f2_yd = {-Two, Two};

Float fn1_f2_epsilon = 1.0e-6f;
Float fn1_f2_delta = 1.0e-20f;
int fn1_f2_nsub = 16;
 
fn1_data_t fn1_f2_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_f2_tag;
    f.descr = fn1_f2_descr;
    f.eval_fp = &fn1_f2_eval_fp;
    f.eval_ia = &fn1_f2_eval_ia;
    f.diff_ia = &fn1_f2_diff_ia;
    f.eval_aa = &fn1_f2_eval_aa;
    f.xd = fn1_f2_xd;
    f.yd = fn1_f2_yd;
    f.epsilon = fn1_f2_epsilon;
    f.delta = fn1_f2_delta;
    f.nsub = fn1_f2_nsub;
    return f;
  }
 
