/* See {mkgr_mark.h} */
/* Last edited on 2020-11-29 18:29:26 by jstolfi */

#define _GNU_SOURCE

#include <vec.h>
#include <affirm.h>

#include <mkgr_mark.h>

mkgr_mark_t mkgr_make_dot(r2_t ctr, double rad, frgb_t color)
  { mkgr_mark_t mk  = (mkgr_mark_t) 
      { .ctr = ctr,
        .cross = FALSE,
        .rad = rad,
        .ang = 0.0,
        .color = color,
        .lwd = 0.0
      };
    return mk;
  }

mkgr_mark_t mkgr_make_circle(r2_t ctr, double rad, double lwd, frgb_t color)
  { demand(lwd > 0.0, "line width must be positive");
    mkgr_mark_t mk  = (mkgr_mark_t) 
      { .ctr = ctr,
        .cross = FALSE,
        .rad = rad,
        .ang = 0.0,
        .color = color,
        .lwd = lwd
      };
    return mk;
  }

mkgr_mark_t mkgr_make_cross(r2_t ctr, double rad, double ang, double lwd, frgb_t color)
  { demand(lwd > 0.0, "line width must be positive");
    mkgr_mark_t mk  = (mkgr_mark_t) 
      { .ctr = ctr,
        .cross = TRUE,
        .rad = rad,
        .ang = ang,
        .color = color,
        .lwd = lwd
      };
    return mk;
  }


vec_typeimpl(mkgr_mark_vec_t, mkgr_mark_vec, mkgr_mark_t);
