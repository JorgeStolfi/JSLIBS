/* See {fn1_gexp.h}. */
/* Last edited on 2013-05-24 04:33:34 by stolfilocal */
 
#include <fn1_gexp.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_gexp_eval_fp (Float x);
Interval fn1_gexp_eval_ia (Interval x);
AAP fn1_gexp_eval_aa (AAP x);
 
/* Last edited on 2001-09-30 02:25:17 by stolfi */

char *fn1_gexp_tag = "gexp";
char *fn1_gexp_descr = "g(x) = exp(x)";

Float fn1_gexp_eval_fp (Float x)
  {
    ROUND_NEAR;
    return (exp(x));
  }

Interval fn1_gexp_eval_ia (Interval x)
  {
    return (ia_exp(x));
  }

AAP fn1_gexp_eval_aa (AAP x)
  {
    /* Interval r; */
    /* r = aa_range(x); */
    return (aa_exp(x));
  }

Interval fn1_gexp_xd = {-Three, One};
Interval fn1_gexp_yd = {-Half, Three+Half};

Float fn1_gexp_epsilon = 1.0e-6;
Float fn1_gexp_delta = 1.0e-20;
int fn1_gexp_nsub = 8;
 
fn1_data_t fn1_gexp_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_gexp_tag;
    f.descr = fn1_gexp_descr;
    f.eval_fp = &fn1_gexp_eval_fp;
    f.eval_ia = &fn1_gexp_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_gexp_eval_aa;
    f.xd = fn1_gexp_xd;
    f.yd = fn1_gexp_yd;
    f.epsilon = fn1_gexp_epsilon;
    f.delta = fn1_gexp_delta;
    f.nsub = fn1_gexp_nsub;
    return f;
  }
 
