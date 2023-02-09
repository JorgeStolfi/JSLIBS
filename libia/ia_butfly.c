/* See ia_butfly.h */
/* Last edited on 2023-02-07 20:42:28 by stolfi */

#include <ia_butfly.h>

#include <flt.h>
#include <ia.h>
#include <ia_trapez.h>
#include <epswr.h>

#include <stdlib.h>
#include <stdio.h>

ia_butfly_t ia_butfly_from_trapez(ia_trapez_t *tp)
  { Interval xe = (Interval){ tp->x.lo, tp->x.lo };
    return (ia_butfly_t) {{ia_trapez_from_box(&xe, &(tp->yxlo)), *tp }};
  }

ia_butfly_t ia_butfly_from_box(Interval *xr, Interval *yr)
  { ia_trapez_t tp = ia_trapez_from_box(xr, yr);
    return ia_butfly_from_trapez(&tp);
  }

ia_butfly_t ia_butfly_from_ia_diff(Interval *xr, Float xm, Interval *yxmr, Interval *dyr)
  { 
    Interval xlo = (Interval){ xr->lo, xm };
    Interval xhi = (Interval){ xm, xr->hi };
    ia_butfly_t bt;
    bt.tp[0] = ia_trapez_from_ia_diff(&xlo, yxmr, dyr, 1);
    bt.tp[1] = ia_trapez_from_ia_diff(&xhi, yxmr, dyr, 0);
    return bt;
  }

void ia_butfly_print(FILE *wr, ia_butfly_t *bt, char *sep)
  { ia_trapez_print(wr, &(bt->tp[0])); 
    if (sep != NULL) { fputs(sep, stderr); }
    ia_trapez_print(wr, &(bt->tp[1])); 
  }

void ia_butfly_draw(epswr_figure_t *eps, Interval *yr, ia_butfly_t *bt)
  { ia_trapez_draw(eps, yr, &(bt->tp[0]));
    ia_trapez_draw(eps, yr, &(bt->tp[1]));
  }
