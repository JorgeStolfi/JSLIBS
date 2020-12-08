/* See {fn1_f3.h}. */
/* Last edited on 2016-12-26 21:26:33 by stolfilocal */
 
#include <fn1_f3.h>
#include <fn1_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
 
Float fn1_f3_eval_fp (Float x);
Interval fn1_f3_eval_ia (Interval x);
Interval fn1_f3_diff_ia (Interval x);
AAP fn1_f3_eval_aa (AAP x);
 
char *fn1_f3_tag = "f3";
char *fn1_f3_descr = 
  "f(x) = sqrt((x-1)^2 + 1) - sqrt((x+1)^2 + 1) - 1/2";

Float fn1_f3_eval_fp (Float x)
  {
    ROUND_NEAR;
    {
      Float xp = x + 1.0f;
      Float xpr = xp * xp + 1.0f;
      Float xps = (Float)sqrt(xpr);
      Float xm = x - 1.0f;
      Float xmr = xm * xm + 1.0f;
      Float xms = (Float)sqrt(xmr);
      Float f = xms - xps - 0.5f;
      return (f);
    }
  }

Interval fn1_f3_eval_ia (Interval x)
  {
    Interval xp = ia_shift(x, +One);
    Interval xpr = ia_shift(ia_mul(xp, xp), One);
    Interval xps = ia_sqrt(xpr);
    Interval xm = ia_shift(x, -One);
    Interval xmr = ia_shift(ia_mul(xm, xm), One);
    Interval xms = ia_sqrt(xmr);
    Interval f = ia_shift(ia_sub(xms, xps), -Half);
    return f;
  }

Interval fn1_f3_diff_ia (Interval x)
  {
    Interval one = {One, One};
    
    Interval xp = ia_shift(x, +One);
    Interval d_xp = one;

    Interval xpr = ia_shift(ia_mul(xp, xp), One);
    Interval d_xpr = ia_mul(ia_scale(xp, Two, One), d_xp);

    Interval xps = ia_sqrt(xpr);
    Interval d_xps = ia_scale(ia_div(d_xpr, xps), One, Two);

    Interval xm = ia_shift(x, -One);
    Interval d_xm = one;

    Interval xmr = ia_shift(ia_mul(xm, xm), One);
    Interval d_xmr = ia_mul(ia_scale(xm, Two, One), d_xm);

    Interval xms = ia_sqrt(xmr);
    Interval d_xms = ia_scale(ia_div(d_xmr, xms), One, Two);
    
    Interval d_f = ia_sub(d_xms, d_xps);

    return d_f;
  }

AAP fn1_f3_eval_aa (AAP x)
  {
    MemP frame = aa_top();
    AAP xp = aa_shift(x, +One);
    AAP xpr = aa_shift(aa_mul(xp, xp), One);
    AAP xps = aa_sqrt(xpr);
    AAP xm = aa_shift(x, -One);
    AAP xmr = aa_shift(aa_mul(xm, xm), One);
    AAP xms = aa_sqrt(xmr);
    AAP f = aa_shift(aa_sub(xms, xps), -Half);
    AAP res = f;
    return (aa_return(frame, res));
  }

Interval fn1_f3_xd = {-Four-Two, Four+Two};
Interval fn1_f3_yd = {-Three, Three};

Float fn1_f3_epsilon = 1.0e-6f;
Float fn1_f3_delta = 1.0e-20f;
int fn1_f3_nsub = 16;
 
fn1_data_t fn1_f3_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_f3_tag;
    f.descr = fn1_f3_descr;
    f.eval_fp = &fn1_f3_eval_fp;
    f.eval_ia = &fn1_f3_eval_ia;
    f.diff_ia = &fn1_f3_diff_ia;
    f.eval_aa = &fn1_f3_eval_aa;
    f.xd = fn1_f3_xd;
    f.yd = fn1_f3_yd;
    f.epsilon = fn1_f3_epsilon;
    f.delta = fn1_f3_delta;
    f.nsub = fn1_f3_nsub;
    return f;
  }
 
