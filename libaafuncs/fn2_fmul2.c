/* See {fn2_fmul2.h}. */
/* Last edited on 2005-09-25 14:23:16 by stolfi */
 
#include <fn2_fmul2.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_fmul2_eval_fp (Float x, Float y);
Interval fn2_fmul2_eval_ia (Interval x, Interval y);
AAP fn2_fmul2_eval_aa (AAP x, AAP y);
 
char *fn2_fmul2_tag = "fmul2";
char *fn2_fmul2_descr = "f = (x - y) * (x + y) - 1";

Float fn2_fmul2_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    { Float u = (x - y);
      Float v = (x + y);
      return (u * v - One);
    }
  }

Interval fn2_fmul2_eval_ia (Interval x, Interval y)
  {
    Interval u = ia_sub(x, y);
    Interval v = ia_add(x, y);
    Interval res = ia_shift(ia_mul(u, v), -One);
    return (res);
  }

AAP fn2_fmul2_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP u = aa_sub(x, y);
    AAP v = aa_add(x, y);
    AAP res = aa_shift(aa_mul(u, v), -One);
    return (aa_return(frame, res));
  }

Interval fn2_fmul2_xd = {-Three, Three};
Interval fn2_fmul2_yd = {-Three, Three};

int fn2_fmul2_fn = 16;
 
fn2_data_t fn2_fmul2_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_fmul2_tag;
    f.descr = fn2_fmul2_descr;
    f.eval_fp = &fn2_fmul2_eval_fp;
    f.eval_ia = &fn2_fmul2_eval_ia;
    f.eval_aa = &fn2_fmul2_eval_aa;
    f.xd = fn2_fmul2_xd;
    f.yd = fn2_fmul2_yd;
    return f;
  }
 
