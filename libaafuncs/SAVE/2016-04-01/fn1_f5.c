/* See {fn1_f5.h}. */
/* Last edited on 2005-09-25 21:03:05 by stolfi */
 
#include <fn1_f5.h>
#include <fn1_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
 
Float fn1_f5_eval_fp (Float x);
Interval fn1_f5_eval_ia (Interval x);
Interval fn1_f5_diff_ia (Interval x);
AAP fn1_f5_eval_aa (AAP x);
 
char *fn1_f5_tag = "f4";
char *fn1_f5_descr = 
  "Q = [+1, -2, +3, -2]; P = [-2, -1, +1, +2];\n" 
  "f(x) = sum_{i=1}^4 Q_i/(x^2 - 2*x*P_i + P_i^2 + 1/2)";

/* Note: we assume that {P[i}*P[i]} has no rounding. */
#define NN 4
static Float Q[NN] = { +1, -2, +3, -2 };
static Float P[NN] = { -2, -2, +1, +2 };

#define Five (5.0)
#define Six (6.0)
#define Fudge Half

Float fn1_f5_eval_fp (Float x)
  { ROUND_NEAR;
    { 
      Float x2 = x*x;
      Float f = 0;
      
      int i;
      for (i = 0; i < NN; i++)
        { Float den = x2 - 2*x*P[i] + P[i]*P[i] + Fudge;
          f = f + Q[i]/den;
        }
      return f;
    }
  }

Interval fn1_f5_eval_ia (Interval x)
  { 
    ROUND_NEAR;
    { 
      Interval x2 = ia_sqr(x);
      Interval f = (Interval){0,0};

      int i;
      for (i = 0; i < NN; i++)
        { Float shf = P[i]*P[i] + Fudge; 
          Interval den = ia_add(x2, ia_shift(ia_scale(x, -2*P[i], One), shf));
          Interval term = ia_scale(ia_inv(den), Q[i], One);
          f = ia_add(f, term);
        }
      return f;
    }
  }

Interval fn1_f5_diff_ia (Interval x)
  { 
    ROUND_NEAR;
    { 
      Interval x2 = ia_sqr(x);
      Interval d_x2 = ia_scale(x, Two, One);

      Interval d_f = (Interval){0,0};

      int i;
      for (i = 0; i < NN; i++)
        { Float shf = P[i]*P[i] + Fudge; 
          Interval den = ia_add(x2, ia_shift(ia_scale(x, -2*P[i], One), shf));
          Interval d_den = ia_shift(d_x2, -2*P[i]);
          Interval d_term = ia_scale(ia_mul(ia_inv(ia_sqr(den)),d_den), -Q[i], One);
          d_f = ia_add(d_f, d_term);
        }
      return d_f;
    }
  }


AAP fn1_f5_eval_aa (AAP x)
  { 
    ROUND_NEAR;
    { 
      MemP frame = aa_top();

      AAP x2 = aa_sqr(x);
      AAP f = aa_zero();

      int i;
      for (i = 0; i < NN; i++)
        { Float shf = P[i]*P[i] + Fudge; 
          AAP den = aa_add(x2, aa_shift(aa_scale(x, -2*P[i], One), shf));
          AAP term = aa_scale(aa_inv(den), Q[i], One);
          f = aa_add(f, term);
        }
      return aa_return(frame, f);
    }
  }

Interval fn1_f5_xd = {-Six, Six};
Interval fn1_f5_yd = {-Six, Six};

Float fn1_f5_epsilon = 5.0e-7;
Float fn1_f5_delta = 1.0e-20;
int fn1_f5_nsub = 16;
 
fn1_data_t fn1_f5_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_f5_tag;
    f.descr = fn1_f5_descr;
    f.eval_fp = &fn1_f5_eval_fp;
    f.eval_ia = &fn1_f5_eval_ia;
    f.diff_ia = &fn1_f5_diff_ia;
    f.eval_aa = &fn1_f5_eval_aa;
    f.xd = fn1_f5_xd;
    f.yd = fn1_f5_yd;
    f.epsilon = fn1_f5_epsilon;
    f.delta = fn1_f5_delta;
    f.nsub = fn1_f5_nsub;
    return f;
  }
 
