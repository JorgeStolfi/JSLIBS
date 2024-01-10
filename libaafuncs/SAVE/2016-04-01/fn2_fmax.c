/* See {fn2_fmax.h}. */
/* Last edited on 2016-04-01 01:14:09 by stolfilocal */
 
#include <fn2_fmax.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_fmax_eval_fp (Float x, Float y);
Interval fn2_fmax_eval_ia (Interval x, Interval y);
AAP fn2_fmax_eval_aa (AAP x, AAP y);
 
char *fn2_fmax_tag = "fmax";
char *fn2_fmax_descr = "f = max(x,y) - 1";

Float fn2_fmax_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    return (float)(FMAX(x, y) - One);
  }

Interval fn2_fmax_eval_ia (Interval x, Interval y)
  {
    Interval res = ia_shift(ia_max(x, y), -One);
    return (res);
  }

AAP fn2_fmax_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP res = aa_shift(aa_max(x, y), -One);
    return (aa_return(frame, res));
  }

Interval fn2_fmax_xd = {-Two, Two};
Interval fn2_fmax_yd = {-Two, Two};

int fn2_fmax_fn = 16;
 
fn2_data_t fn2_fmax_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_fmax_tag;
    f.descr = fn2_fmax_descr;
    f.eval_fp = &fn2_fmax_eval_fp;
    f.eval_ia = &fn2_fmax_eval_ia;
    f.eval_aa = &fn2_fmax_eval_aa;
    f.xd = fn2_fmax_xd;
    f.yd = fn2_fmax_yd;
    return f;
  }
 
