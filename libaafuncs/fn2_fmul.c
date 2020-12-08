/* See {fn2_fmul.h}. */
/* Last edited on 2005-09-25 14:23:10 by stolfi */
 
#include <fn2_fmul.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_fmul_eval_fp (Float x, Float y);
Interval fn2_fmul_eval_ia (Interval x, Interval y);
AAP fn2_fmul_eval_aa (AAP x, AAP y);
 
char *fn2_fmul_tag = "fmul";
char *fn2_fmul_descr = "f = x * y - 1";

Float fn2_fmul_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    return (x * y - One);
  }

Interval fn2_fmul_eval_ia (Interval x, Interval y)
  {
    Interval res = ia_shift(ia_mul(x, y), -One);
    return (res);
  }

AAP fn2_fmul_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP res = aa_shift(aa_mul(x, y), -One);
    return (aa_return(frame, res));
  }

Interval fn2_fmul_xd = {-Two, Two};
Interval fn2_fmul_yd = {-Two, Two};

int fn2_fmul_fn = 16;
 
fn2_data_t fn2_fmul_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_fmul_tag;
    f.descr = fn2_fmul_descr;
    f.eval_fp = &fn2_fmul_eval_fp;
    f.eval_ia = &fn2_fmul_eval_ia;
    f.eval_aa = &fn2_fmul_eval_aa;
    f.xd = fn2_fmul_xd;
    f.yd = fn2_fmul_yd;
    return f;
  }
 
