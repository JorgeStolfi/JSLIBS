/* See {fn2_f8.h}. */
/* Last edited on 2016-04-01 01:13:38 by stolfilocal */
 
#include <fn2_f8.h>
#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
 
Float fn2_f8_eval_fp (Float x, Float y);
Interval fn2_f8_eval_ia (Interval x, Interval y);
AAP fn2_f8_eval_aa (AAP x, AAP y);
 
char *fn2_f8_tag = "f8";
char *fn2_f8_descr = 
  "f(x,y) = sum_i q_i / dist((x,y), (x_i, y_i))";

#define x1 (+0.500f)
#define y1 (+0.500f)
#define q1 (+1.000f)

#define x2 (-0.500f)
#define y2 (+0.250f)
#define q2 (-0.875f)

#define x3 (+0.500f)
#define y3 (-0.125f)
#define q3 (-0.750f)

#define x4 (-0.500f)
#define y4 (-0.375f)
#define q4 (+0.625f)

Float fn2_f8_d_fp (Float x, Float y, Float x0, Float y0);
Interval fn2_f8_d_ia (Interval x, Interval y, Float x0, Float y0);
AAP fn2_f8_d_aa (AAP x, AAP y, Float x0, Float y0);

Float fn2_f8_d_fp (Float x, Float y, Float x0, Float y0)
  { 
    ROUND_NEAR;
    {
      Float dx = x - x0;
      Float dy = y - y0;
      Float d = sqrt(dx*dx + dy*dy);
      return(d);
    }
  }
  
Float fn2_f8_eval_fp (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float p1 = q1 / fn2_f8_d_fp (x, y, x1, y1);
      Float p2 = q2 / fn2_f8_d_fp (x, y, x2, y2);
      Float p3 = q3 / fn2_f8_d_fp (x, y, x3, y3);
      Float p4 = q4 / fn2_f8_d_fp (x, y, x4, y4);
      Float res = p1 + p2 + p3 + p4;
      return (res);
    }
  }

Interval fn2_f8_d_ia (Interval x, Interval y, Float x0, Float y0)
  { 
    Interval dx = ia_shift(x, -x0);
    Interval dy = ia_shift(y, -y0);
    Interval d = ia_sqrt(ia_add(ia_sqr(dx), ia_sqr(dy)));
    return(d);
  }
  
Interval fn2_f8_eval_ia (Interval x, Interval y)
  {
    Interval p1 = ia_div(ia_const(q1, Zero), fn2_f8_d_ia (x, y, x1, y1));
    Interval p2 = ia_div(ia_const(q2, Zero), fn2_f8_d_ia (x, y, x2, y2));
    Interval p3 = ia_div(ia_const(q3, Zero), fn2_f8_d_ia (x, y, x3, y3));
    Interval p4 = ia_div(ia_const(q4, Zero), fn2_f8_d_ia (x, y, x4, y4));
    Interval res = ia_add(ia_add(p1, p2), ia_add(p3, p4));
    return (res);
  }

AAP fn2_f8_d_aa (AAP x, AAP y, Float x0, Float y0)
  { 
    AAP dx = aa_shift(x, -x0);
    AAP dy = aa_shift(y, -y0);
    AAP d = aa_sqrt(aa_add(aa_sqr(dx), aa_sqr(dy)));
    return(d);
  }
  
AAP fn2_f8_eval_aa (AAP x, AAP y)
  {
    MemP frame = aa_top();
    AAP p1 = aa_div(aa_const(q1, Zero), fn2_f8_d_aa (x, y, x1, y1));
    AAP p2 = aa_div(aa_const(q2, Zero), fn2_f8_d_aa (x, y, x2, y2));
    AAP p3 = aa_div(aa_const(q3, Zero), fn2_f8_d_aa (x, y, x3, y3));
    AAP p4 = aa_div(aa_const(q4, Zero), fn2_f8_d_aa (x, y, x4, y4));
    AAP res = aa_add(aa_add(p1, p2), aa_add(p3, p4));
    return (aa_return(frame, res));
  }

Interval fn2_f8_xd = {-Two, Two};
Interval fn2_f8_yd = {-Two, Two};

int fn2_f8_fn = 64;
 
fn2_data_t fn2_f8_get_data(void)
  { 
    fn2_data_t f;
    f.tag = fn2_f8_tag;
    f.descr = fn2_f8_descr;
    f.eval_fp = &fn2_f8_eval_fp;
    f.eval_ia = &fn2_f8_eval_ia;
    f.eval_aa = &fn2_f8_eval_aa;
    f.xd = fn2_f8_xd;
    f.yd = fn2_f8_yd;
    return f;
  }
 
