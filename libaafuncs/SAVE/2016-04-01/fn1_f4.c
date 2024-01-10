/* See {fn1_f4.h}. */
/* Last edited on 2005-09-25 19:12:26 by stolfi */
 
#include <fn1_f4.h>
#include <fn1_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
 
Float fn1_f4_eval_fp (Float x);
Interval fn1_f4_eval_ia (Interval x);
Interval fn1_f4_diff_ia (Interval x);
AAP fn1_f4_eval_aa (AAP x);
 
char *fn1_f4_tag = "f4";
char *fn1_f4_descr = 
  "Q = [+1, -2, +3, -2]; P = [-2, -1, +1, +2];\n" 
  "f(x) = sum_{i=1}^4 Q_i/((x - P_i)^2 + 1/2)";

#define Qa (+One)
#define Qb (-Two)
#define Qc (+Three)
#define Qd (-Two)
#define Five (5.0)
#define Six (6.0)
#define Fudge Half

Float fn1_f4_eval_fp (Float x)
  { ROUND_NEAR;
    { Float xa = (x + 2.0)*(x + 2.0) + Fudge;
      Float xb = (x + 1.0)*(x + 1.0) + Fudge;
      Float xc = (x - 1.0)*(x - 1.0) + Fudge;
      Float xd = (x - 2.0)*(x - 2.0) + Fudge;
      
      Float fa = Qa/xa;
      Float fb = Qb/xb;
      Float fc = Qc/xc;
      Float fd = Qd/xd;
      
      Float f = fa + fb + fc + fd;
      return (f);
    }
  }

Interval fn1_f4_eval_ia (Interval x)
  { 
    Interval xa = ia_shift(ia_sqr(ia_shift(x, +Two)), Fudge);
    Interval xb = ia_shift(ia_sqr(ia_shift(x, +One)), Fudge);
    Interval xc = ia_shift(ia_sqr(ia_shift(x, -One)), Fudge);
    Interval xd = ia_shift(ia_sqr(ia_shift(x, -Two)), Fudge);
    Interval fa = ia_scale(ia_inv(xa), Qa, One);
    Interval fb = ia_scale(ia_inv(xb), Qb, One);
    Interval fc = ia_scale(ia_inv(xc), Qc, One);
    Interval fd = ia_scale(ia_inv(xd), Qd, One);
    Interval f = ia_add(ia_add(fa, fb), ia_add(fc, fd));
    return f;
  }

Interval fn1_f4_diff_ia (Interval x)
  { 
    Interval xa = ia_shift(ia_sqr(ia_shift(x, +Two)), Fudge);
    Interval d_xa = ia_scale(ia_shift(x, +Two), Two, One);

    Interval xb = ia_shift(ia_sqr(ia_shift(x, +One)), Fudge);
    Interval d_xb = ia_scale(ia_shift(x, +One), Two, One);

    Interval xc = ia_shift(ia_sqr(ia_shift(x, -One)), Fudge);
    Interval d_xc = ia_scale(ia_shift(x, -One), Two, One);

    Interval xd = ia_shift(ia_sqr(ia_shift(x, -Two)), Fudge);
    Interval d_xd = ia_scale(ia_shift(x, -Two), Two, One);

    Interval d_fa = ia_scale(ia_div(d_xa, ia_sqr(xa)), Qa, -One);

    Interval d_fb = ia_scale(ia_div(d_xb, ia_sqr(xb)), Qb, -One);

    Interval d_fc = ia_scale(ia_div(d_xc, ia_sqr(xc)), Qc, -One);

    Interval d_fd = ia_scale(ia_div(d_xd, ia_sqr(xd)), Qd, -One);

    Interval d_f = ia_add(ia_add(d_fa, d_fb), ia_add(d_fc, d_fd));

    return d_f;
  }


AAP fn1_f4_eval_aa (AAP x)
  { ROUND_NEAR;
    { MemP frame = aa_top();

      AAP xa = aa_shift(aa_sqr(aa_shift(x, +Two)), Fudge);
      AAP xb = aa_shift(aa_sqr(aa_shift(x, +One)), Fudge);
      AAP xc = aa_shift(aa_sqr(aa_shift(x, -One)), Fudge);
      AAP xd = aa_shift(aa_sqr(aa_shift(x, -Two)), Fudge);

      AAP fa = aa_scale(aa_inv(xa), Qa, One);
      AAP fb = aa_scale(aa_inv(xb), Qb, One);
      AAP fc = aa_scale(aa_inv(xc), Qc, One);
      AAP fd = aa_scale(aa_inv(xd), Qd, One);

      AAP f = aa_add(aa_add(fa, fb), aa_add(fc, fd));

      AAP res = f;
      return (aa_return(frame, res));
    }
  }

Interval fn1_f4_xd = {-Six, Six};
Interval fn1_f4_yd = {-Six, Six};

Float fn1_f4_epsilon = 5.0e-7;
Float fn1_f4_delta = 1.0e-20;
int fn1_f4_nsub = 16;
 
fn1_data_t fn1_f4_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_f4_tag;
    f.descr = fn1_f4_descr;
    f.eval_fp = &fn1_f4_eval_fp;
    f.eval_ia = &fn1_f4_eval_ia;
    f.diff_ia = &fn1_f4_diff_ia;
    f.eval_aa = &fn1_f4_eval_aa;
    f.xd = fn1_f4_xd;
    f.yd = fn1_f4_yd;
    f.epsilon = fn1_f4_epsilon;
    f.delta = fn1_f4_delta;
    f.nsub = fn1_f4_nsub;
    return f;
  }
 
