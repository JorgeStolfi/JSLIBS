/* See frgb_interp_vis.h */
/* Last edited on 2023-03-07 13:57:54 by stolfi */

/* Copyright (C) 2023 by Jorge Stolfi, the University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <frgb.h>
#include <frgb_ops.h>

#include <frgb_interp_vis.h>

frgb_t frgb_interp_vis(double r, frgb_t *p0, frgb_t *p1, bool_t logY)
  { frgb_t q0 = *p0; frgb_to_HTY(&q0);
    frgb_t q1 = *p1; frgb_to_HTY(&q1);
    frgb_t p = frgb_interp_vis_HTY(r, &q0, &q1, logY);
    frgb_from_HTY(&p);
    frgb_clip_rgb(&p);
    return p;
  }
  
frgb_t frgb_interp_vis_HTY(double r, frgb_t *q0, frgb_t *q1, bool_t logY)
  { demand((r >= 0) && (r <= 1), "fraction {r} should be in {[0_1]}");

    double H0 = q0->c[0], T0 = q0->c[1], Y0 = q0->c[2]; 
    demand((Y0 >= 0) && (T0 >= 0), "negative {Y} or {T} in {q0}");
    double H1 = q1->c[0], T1 = q1->c[1], Y1 = q1->c[2]; 
    if ((Y1 < 0) || (T1 < 0)) { fprintf(stderr, "Y1 = %24.16e T1 = %24.26e\n", Y1, T1); }
    demand((Y1 >= 0) && (T1 >= 0), "negative {Y} or {T} in {q1}");

    /* If either color is gray, use the other color's hue: */
    if (T0 <= 0) { H0 = H1; }
    if (T1 <= 0) { H1 = H0; }

    double Y;
    if (logY)
      { /* Interpolate the brightness in biased log scale: */
        double Ybias = 0.010;
        Y = exp((1-r)*log(Y0 + Ybias) + r*log(Y1 + Ybias)) - Ybias;
      }
    else
      { /* Interpolate the brightness affinely: */
        Y = (1-r)*Y0 + r*Y1;
      }
    Y = fmin(1, fmax(0, Y));

    frgb_t q;
    if (Y <= 0)
      { q = (frgb_t){{ 0, 0, 0 }}; }
    else if (Y >= 1)
      { q = (frgb_t){{ 1, 1, 1 }}; }
    else
      { /* Interpolate angular {UV} hues {H0,H1} in linear scale */
        double H = (1-r)*H0 + r*H1;

        /* Interpolate relative saturations {Y0,T1} in biased log scale: */
        double Tbias = 0.010;
        double T = exp((1-r)*log(T0 + Tbias) + r*log(T1 + Tbias)) - Tbias;
        T = fmax(0, T);

        /* Repackage:  */
        q = (frgb_t){{ (float)H, (float)T, (float)Y }};
      }

    return q;
  }
