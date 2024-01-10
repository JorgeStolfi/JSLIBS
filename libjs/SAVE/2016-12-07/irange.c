/* See irange.h */
/* Last edited on 2007-05-09 08:18:07 by stolfi */

#include <irange.h>
#include <stdint.h>
#include <stdlib.h>

irange_t irange_join(irange_t *u, irange_t *v)
  { int32_t ulo = ILO(*u), uhi = IHI(*u);
    int32_t vlo = ILO(*v), vhi = IHI(*v);
    irange_t w;
    ILO(w) = (ulo < vlo ? ulo : vlo);
    IHI(w) = (uhi > vhi ? uhi : vhi);
    return w;
  }

irange_t irange_meet(irange_t *u, irange_t *v)
  { int32_t ulo = ILO(*u), uhi = IHI(*u);
    int32_t vlo = ILO(*v), vhi = IHI(*v);
    irange_t w;
    ILO(w) = (ulo > vlo ? ulo : vlo);
    IHI(w) = (uhi < vhi ? uhi : vhi);
    return w;
  }

void irange_widen(irange_t *r, int32_t margin)
  { ILO(*r) = ILO(*r) - margin;
    IHI(*r) = IHI(*r) + margin;
  }
