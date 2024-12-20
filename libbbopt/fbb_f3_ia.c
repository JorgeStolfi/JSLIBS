/* Defines a target function F for bbopt2.c */
/* Last edited on 2024-12-05 10:21:41 by stolfi */
  
#include <stdint.h>

#include <ia.h>

#include <bbgoal.h>
#include <fbb_f3_ia.h>

/* PROTOTYPES */

int32_t fbb_f3_ia_dim = 2;
char *fbb_f3_ia_tag = "f3_ia";
char *fbb_f3_ia_descr = "(IA) u = x - 2/7; v = y - 3/7; u^2 + v^2 + u*v + 1";

Float fbb_f3_ia_eval_fp(Float *x);
Interval fbb_f3_ia_eval_ia(Interval *xr);
void fbb_f3_ia_true_opt(Interval *xr, Interval *sr);

/* IMPLEMENTATIONS */

#define fbb_f3_XMINN (2.0)
#define fbb_f3_XMIND (7.0)
#define fbb_f3_XMIN (fbb_f3_XMINN/fbb_f3_XMIND)

#define fbb_f3_YMINN (3.0)
#define fbb_f3_YMIND (7.0)
#define fbb_f3_YMIN (fbb_f3_YMINN/fbb_f3_YMIND)

Float fbb_f3_ia_eval_fp(Float *x)
  { double dx = x[0] - fbb_f3_XMIN;
    double dy = x[0] - fbb_f3_YMIN;
    double xx = dx*dx;
    double xy = dx*dy;
    double yy = dy*dy;
    return (Float)(xx+xy+yy + 1.0);
  }

Interval fbb_f3_ia_eval_ia(Interval *xr)
  { Interval x = xr[0];
    Interval y = xr[1];
    
    Interval dx = ia_shift(x, (Float)(-2.0/7.0));
    Interval dy = ia_shift(y, (Float)(-3.0/7.0));

    Interval xx = ia_sqr(dx);
    Interval yy = ia_sqr(dy);
    Interval xy = ia_mul(dx, dy);

    Interval f = ia_shift(ia_add(ia_add(xx, yy), xy), +1.0);
    
    return f;
  }

void fbb_f3_ia_true_opt(Interval *xr, Interval *sr)
  { /* Should take `xr' into account... */
    Interval mr[2];
    mr[0] = ia_scale(ia_const(One, Zero), fbb_f3_XMINN, fbb_f3_XMIND);
    mr[1] = ia_scale(ia_const(One, Zero), fbb_f3_YMINN, fbb_f3_YMIND);
    int32_t i;
    for (i = 0; i < fbb_f3_ia_dim; i++)
      { Float mlo = mr[i].lo, mhi = mr[i].hi;
        Float xlo = xr[i].lo, xhi = xr[i].hi;
        Float slo, shi;
        if (mhi <= xlo)
          { slo = shi = xlo; }
        else if (mlo >= xhi)
          { slo = shi = xhi; }
        else
          { slo = (mlo < xlo ? xlo : mlo);
            shi = (mhi > xhi ? xhi : mhi);
          }
        sr[i] = (Interval) { slo, shi };
      }
  }

bbgoal_data_t fbb_f3_ia_get_data(void)
  {
    bbgoal_data_t f;
    f.dim = fbb_f3_ia_dim;
    f.eval_fp = &fbb_f3_ia_eval_fp;
    f.eval_ia = &fbb_f3_ia_eval_ia;
    f.true_opt = &fbb_f3_ia_true_opt;
    f.plot_range = NULL;
    f.tag = fbb_f3_ia_tag;
    f.descr = fbb_f3_ia_descr;
    return f;
  }

