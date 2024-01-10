/* See {fn1_gdiv.h}. */
/* Last edited on 2005-09-25 15:47:31 by stolfi */
 
#include <fn1_gdiv.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_gdiv_eval_fp (Float x);
Interval fn1_gdiv_eval_ia (Interval x);
AAP fn1_gdiv_eval_aa (AAP x);
 
char *fn1_gdiv_tag = "gdiv";
char *fn1_gdiv_descr = 
  "u = sqrt(x^2 + 16);\n"
  "v = sqrt(x^2 + 1);\n"
  "g = u/v";

Float fn1_gdiv_eval_fp (Float x)
  {
    ROUND_NEAR;
    return (sqrt(x*x+Four*Four)/sqrt(x*x+One));
  }

Interval fn1_gdiv_eval_ia (Interval x)
  {
    Interval u = ia_sqrt(ia_shift(ia_sqr(x), Four*Four));
    Interval v = ia_sqrt(ia_shift(ia_sqr(x), One));
    Interval res = ia_div(u,v);
    return (res);
  }

AAP fn1_gdiv_eval_aa (AAP x)
  {
    AAP u = aa_sqrt(aa_shift(aa_sqr(x), Four*Four));
    AAP v = aa_sqrt(aa_shift(aa_sqr(x), One));
    AAP res = aa_div(u,v);
    return (res);
  }

Interval fn1_gdiv_xd = {-One, Two*Four - One};
Interval fn1_gdiv_yd = {-Two, Two*Four};

Float fn1_gdiv_epsilon = 1.0e-6;
Float fn1_gdiv_delta = 1.0e-20;
int fn1_gdiv_nsub = 12;
 
fn1_data_t fn1_gdiv_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_gdiv_tag;
    f.descr = fn1_gdiv_descr;
    f.eval_fp = &fn1_gdiv_eval_fp;
    f.eval_ia = &fn1_gdiv_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_gdiv_eval_aa;
    f.xd = fn1_gdiv_xd;
    f.yd = fn1_gdiv_yd;
    f.epsilon = fn1_gdiv_epsilon;
    f.delta = fn1_gdiv_delta;
    f.nsub = fn1_gdiv_nsub;
    return f;
  }
 
