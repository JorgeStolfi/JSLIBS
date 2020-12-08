/* See {fn1_f1.h}. */
/* Last edited on 2016-12-26 21:34:41 by stolfilocal */
 
#include <fn1_f1.h>
#include <fn1_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
 
Float fn1_f1_eval_fp (Float x);
Interval fn1_f1_eval_ia (Interval x);
Interval fn1_f1_diff_ia (Interval x);
AAP fn1_f1_eval_aa (AAP x);

char *fn1_f1_tag = "f1";
char *fn1_f1_descr = 
  "f(x) = sqrt(x^2 - x + 1/2)/sqrt(x^2 + 1/2) - 1";
 
Float fn1_f1_eval_fp (Float x)
  { ROUND_NEAR;
    return (Float)(sqrt(x*x - x + 0.50)/sqrt(x*x + 0.50) - 1.0);
  }

Interval fn1_f1_eval_ia (Interval x)
  { 
    Interval hlf = {Half, Half};
    Interval x2  = ia_sqr(x);
    Interval dif = ia_sub(x2, x);
    Interval sum1 = ia_add(dif, hlf);
    Interval sum2 = ia_add(x2, hlf);
    Interval sqrt1 = ia_sqrt(sum1);
    Interval sqrt2 = ia_sqrt(sum2);
    Interval quo = ia_div(sqrt1, sqrt2);
    Interval f = ia_shift(quo, -One);
    return f;
  }

Interval fn1_f1_diff_ia (Interval x)
  { 
    Interval hlf = {Half, Half};

    Interval x2  = ia_sqr(x);
    Interval d_x2 = ia_scale(x, Two, One);

    Interval dif = ia_sub(x2, x);
    Interval d_dif = ia_shift(d_x2, -One);

    Interval sum1 = ia_add(dif, hlf);
    Interval d_sum1 = d_dif;

    Interval sum2 = ia_add(x2, hlf);
    Interval d_sum2 = d_x2;
    
    Interval sqrt1 = ia_sqrt(sum1);
    Interval d_sqrt1 = ia_mul(ia_scale(ia_inv(sqrt1), Half, One), d_sum1);
    
    Interval sqrt2 = ia_sqrt(sum2);
    Interval d_sqrt2 = ia_mul(ia_scale(ia_inv(sqrt2), Half, One), d_sum2);
    
    Interval quo = ia_div(sqrt1, sqrt2);
    Interval d_quo = ia_meet(
      ia_div(ia_sub(d_sqrt1, ia_mul(quo, d_sqrt2)), sqrt2), 
      ia_div(ia_sub(ia_mul(d_sqrt1, sqrt2), ia_mul(sqrt1, d_sqrt2)), ia_sqr(sqrt2))
    );
    
    Interval d_f = d_quo;

    return d_f;
  }

AAP fn1_f1_eval_aa (AAP x)
  { MemP frame = aa_top();
    AAP hlf = aa_const(Half, Zero);
    AAP x2  = aa_sqr(x);
    AAP dif = aa_sub(x2, x);
    AAP sum1 = aa_add(dif, hlf);
    AAP sum2 = aa_add(x2, hlf);
    AAP quo = aa_div(aa_sqrt(sum1), aa_sqrt(sum2));
    AAP res = aa_shift(quo, -One);
    return (aa_return(frame, res));
  }

Interval fn1_f1_xd = {-Two, Two};
Interval fn1_f1_yd = {-(One+Half), One+Half};

Float fn1_f1_epsilon = 1.0e-6f;
Float fn1_f1_delta = 1.0e-6f;
int fn1_f1_nsub = 16;
 
fn1_data_t fn1_f1_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_f1_tag;
    f.descr = fn1_f1_descr;
    f.eval_fp = &fn1_f1_eval_fp;
    f.eval_ia = &fn1_f1_eval_ia;
    f.diff_ia = &fn1_f1_diff_ia;
    f.eval_aa = &fn1_f1_eval_aa;
    f.xd = fn1_f1_xd;
    f.yd = fn1_f1_yd;
    f.epsilon = fn1_f1_epsilon;
    f.delta = fn1_f1_delta;
    f.nsub = fn1_f1_nsub;
    return f;
  }
 
