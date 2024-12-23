/* See pst_geom.h */
/* Last edited on 2024-12-22 12:04:21 by stolfi */

#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <r2.h>
#include <r3.h> 
#include <affirm.h>

#include <pst_geom.h>

#define Pr fprintf
#define Er stderr

#define X c[0]
#define Y c[1]
#define Z c[2]
  /* Cartesian coordinates of {r2_t} or {r3_t}. */

/* UTILITIES */

void pst_geom_clip_dir(r3_t *udir, r3_t *sdir, double ard)      
  { double ucos = r3_dot(sdir, udir);
    if (ucos < cos(ard))
      { r3_t para, perp; /* Comps. of {udir} parallel and perpendic. to {sdir}. */
        r3_decomp(udir, sdir, &para, &perp);
        while (r3_L_inf_norm(&perp) == 0.0)
          { /* Set {perp} to a random direction perpendicular to {sdir}. */
            r3_throw_ball(&perp);
            r3_decomp(&perp, sdir, &para, &perp);
          }
        /* Normalize {perp} to unit length. */
        (void)r3_dir(&perp, &perp);
        /* Set {udir} to point {ard} radians away from {sdir}: */
        double cr = cos(ard), sr = sin(ard);
        r3_mix(cr, sdir, sr, &perp, udir);
      }
  }

