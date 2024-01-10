/* See {fn2_fadd.h}. */
/* Last edited on 2005-09-25 14:22:49 by stolfi */
 
#include <fn2_fadd.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_fadd_eval_fp (Float x, Float y);
Interval fn2_fadd_eval_ia (Interval x, Interval y);
AAP fn2_fadd_eval_aa (AAP x, AAP y);
 
char *fn2_fadd_tag = "fadd";
char *fn2_fadd_descr = "f = x + y - 1";

Float fn2_fadd_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    return (x + y - One);
  }

Interval fn2_fadd_eval_ia (Interval x, Interval y)
  {
    Interval res = ia_shift(ia_add(x, y), -One);
    return (res);
  }

AAP fn2_fadd_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP res = aa_shift(aa_add(x, y), -One);
    return (aa_return(frame, res));
  }

Interval fn2_fadd_xd = {-Two, Two};
Interval fn2_fadd_yd = {-Two, Two};

int fn2_fadd_fn = 16;
 
fn2_data_t fn2_fadd_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_fadd_tag;
    f.descr = fn2_fadd_descr;
    f.eval_fp = &fn2_fadd_eval_fp;
    f.eval_ia = &fn2_fadd_eval_ia;
    f.eval_aa = &fn2_fadd_eval_aa;
    f.xd = fn2_fadd_xd;
    f.yd = fn2_fadd_yd;
    return f;
  }
 
