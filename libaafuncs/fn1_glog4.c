/* See {fn1_glog4.h}. */
/* Last edited on 2016-12-26 21:34:08 by stolfilocal */
 
#include <fn1_glog4.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn1_functions.h>
 
Float fn1_glog4_eval_fp (Float x);
Interval fn1_glog4_eval_ia (Interval x);
AAP fn1_glog4_eval_aa (AAP x);
 
char *fn1_glog4_tag = "glog4";
char *fn1_glog4_descr = 
  "u = x - 1;\n"
  "u2 = y^2;\n"
  "u3 = u*u2;\n"
  "u4 = u2^2;\n"
  "g(x) = u - u2/2 + u3/3 - u4/4";

#define Six (6.0)
#define Twenty (20.0)

Float fn1_glog4_eval_fp (Float x)
  {
    ROUND_NEAR;
    { 
      Float u = x - One;
      Float u2 = u*u;
      Float u3 = u*u2;
      Float u4 = u2*u2;
      Float t1 = u;
      Float t2 = u2/Two;
      Float t3 = u3/Three;
      Float t4 = u4/Four;
      Float res = t1 - t2 + t3 - t4;
      return (res);
    }
  }

Interval fn1_glog4_eval_ia (Interval x)
  {
    Interval u = ia_shift(x, -One);
    Interval u2 = ia_sqr(u);
    Interval u3 = ia_mul(u, u2);
    Interval u4 = ia_sqr(u2);
    Interval t1 = u;
    Interval t2 = ia_scale(u2, One, Two);
    Interval t3 = ia_scale(u3, One, Three);
    Interval t4 = ia_scale(u4, One, Four);
    Interval res = ia_add(ia_sub(t1, t2), ia_sub(t3, t4));
    return (res);
  }

AAP fn1_glog4_eval_aa (AAP x)
  { 
    MemP frame = aa_top();
    AAP u = aa_shift(x, -One);
    AAP u2 = aa_sqr(u);
    AAP u3 = aa_mul(u, u2);
    AAP u4 = aa_sqr(u2);
    AAP t1 = u;
    AAP t2 = aa_scale(u2, One, Two);
    AAP t3 = aa_scale(u3, One, Three);
    AAP t4 = aa_scale(u4, One, Four);
    AAP res = aa_add(aa_sub(t1, t2), aa_sub(t3, t4));
    return (aa_return(frame, res));
  }

Interval fn1_glog4_xd = {-One, Three};
Interval fn1_glog4_yd = {-Three, One};

Float fn1_glog4_epsilon = 1.0e-6f;
Float fn1_glog4_delta = 1.0e-20f;
int fn1_glog4_nsub = 16;
 
fn1_data_t fn1_glog4_get_data(void)
  { 
    fn1_data_t f;
    f.tag = fn1_glog4_tag;
    f.descr = fn1_glog4_descr;
    f.eval_fp = &fn1_glog4_eval_fp;
    f.eval_ia = &fn1_glog4_eval_ia;
    f.diff_ia = NULL;
    f.eval_aa = &fn1_glog4_eval_aa;
    f.xd = fn1_glog4_xd;
    f.yd = fn1_glog4_yd;
    f.epsilon = fn1_glog4_epsilon;
    f.delta = fn1_glog4_delta;
    f.nsub = fn1_glog4_nsub;
    return f;
  }
 
