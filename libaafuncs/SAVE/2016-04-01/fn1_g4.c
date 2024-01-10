/* See {fn1_g4.h}. */
/* Last edited on 2005-09-25 15:46:27 by stolfi */
 
#include <fn1_g4.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_g4_eval_fp (Float x);
Interval fn1_g4_eval_ia (Interval x);
AAP fn1_g4_eval_aa (AAP x);
 
char *fn1_g4_tag = "g4";
char *fn1_g4_descr = 
  "u = x - 1; v = x + 1;\n"
  "a = u^2 + 1/2; b = v2 + 1/2;\n"
  "g = 1/a + 1/b";

Float fn1_g4_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float u = x - One;
      Float v = x + One;
      Float u2 = u*u;
      Float v2 = v*v;
      Float a = u2 + Half;
      Float b = v2 + Half;
      Float res = 1/a + 1/b;
      return (res);
    }
  }

Interval fn1_g4_eval_ia (Interval x)
  {
    Interval u2 = ia_sqr(ia_shift(x, -One));
    Interval v2 = ia_sqr(ia_shift(x, +One));
    Interval a = ia_shift(u2, Half);
    Interval b = ia_shift(v2, Half);
    Interval res = ia_add(ia_inv(a), ia_inv(b));
    return (res);
  }

AAP fn1_g4_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP u2 = aa_sqr(aa_shift(x, -One));
    AAP v2 = aa_sqr(aa_shift(x, +One));
    AAP a = aa_shift(u2, Half);
    AAP b = aa_shift(v2, Half);
    AAP res = aa_add(aa_inv(a), aa_inv(b));
    return (aa_return(frame, res));
  }

Interval fn1_g4_xd = {-Three, Three};
Interval fn1_g4_yd = {-Three, Three};

Float fn1_g4_epsilon = 1.0e-6;
Float fn1_g4_delta = 1.0e-20;
int fn1_g4_nsub = 32;
 
fn1_data_t fn1_g4_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_g4_tag;
    f.descr = fn1_g4_descr;
    f.eval_fp = &fn1_g4_eval_fp;
    f.eval_ia = &fn1_g4_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_g4_eval_aa;
    f.xd = fn1_g4_xd;
    f.yd = fn1_g4_yd;
    f.epsilon = fn1_g4_epsilon;
    f.delta = fn1_g4_delta;
    f.nsub = fn1_g4_nsub;
    return f;
  }
 
