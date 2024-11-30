/* See fn2_zf_flt.h */
/* Last edited on 2023-02-18 04:28:52 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>

#include <flt.h>
#include <affirm.h>
#include <ia.h>
#include <epswr.h>

#include <fn2_zf_flt.h>

/*** PROTOTYPES FOR INTERNAL ROUTINES ***/

void fn2_zf_flt_in_triangle(
    epswr_figure_t *eps,
    double xa, double ya, Float za,
    double xb, double yb, Float zb,
    double xc, double yc, Float zc
  );
  /* Draws intersection of triangle abc with plane z=0. */
  /* Includes edges ab and bc but not ac.               */

/*** IMPLEMENTATIONS ***/

void fn2_zf_flt_plot(
    epswr_figure_t *eps,
    Float f (Float x, Float y),
    Interval xd,
    Interval yd,
    int32_t m
  )
  {
    Interval xv, yv;
    Float xm, ym;
    Float f00, f01, f10, f11, fmm;

    ROUND_NEAR;
    epswr_comment(eps, "Plot of actual zeros");

    epswr_set_pen(eps, 0.5, 0.0, 0.0, 0.30, 0.0, 0.0);

    for (uint32_t yi = 0;  yi<m; yi++)
      { for (uint32_t xi = 0;  xi<m; xi++)
          { ROUND_DOWN;
            xv.lo = xd.lo + ((xd.hi - xd.lo)*(Float)xi)/(Float)m;
            yv.lo = yd.lo + ((yd.hi - yd.lo)*(Float)yi)/(Float)m;

            ROUND_UP;
            xv.hi = xd.lo + ((xd.hi - xd.lo)*(Float)(xi+1))/(Float)m;
            yv.hi = yd.lo + ((yd.hi - yd.lo)*(Float)(yi+1))/(Float)m;

            ROUND_NEAR;
            xm = (xv.lo + xv.hi)/Two;
            ym = (yv.lo + yv.hi)/Two;

            f00 = f(xv.lo, yv.lo);
            f01 = f(xv.lo, yv.hi);
            f10 = f(xv.hi, yv.lo);
            f11 = f(xv.hi, yv.hi);
            fmm = f(xm, ym);

            fn2_zf_flt_in_triangle(eps,
              xv.lo,  yv.lo,  f00,
              xv.hi,  yv.lo,  f10,
              xm,     ym,     fmm
            );

            fn2_zf_flt_in_triangle(eps,
              xv.hi,  yv.lo,  f10,
              xv.hi,  yv.hi,  f11,
              xm,     ym,     fmm
            );

            fn2_zf_flt_in_triangle(eps,
              xv.hi,  yv.hi,  f11,
              xv.lo,  yv.hi,  f01,
              xm,     ym,     fmm
            );

            fn2_zf_flt_in_triangle(eps,
              xv.lo,  yv.hi,  f01,
              xv.lo,  yv.lo,  f00,
              xm,     ym,     fmm
            );

          }
      }
  }

void fn2_zf_flt_in_triangle(
    epswr_figure_t *eps,
    double xa, double ya, Float za,
    double xb, double yb, Float zb,
    double xc, double yc, Float zc
  )
  { if (za == Zero && zb == Zero && zc == Zero)
      {
        epswr_set_fill_color(eps, 0.0, 0.0, 0.0);
        epswr_triangle(eps, xa, ya, xb, yb, xc, yc, TRUE, FALSE);
      }
    if (za == Zero && zb == Zero)
      {
        epswr_segment(eps, xa, ya, xb, yb);
      }
    if (zb == Zero && zc == Zero)
      {
        epswr_segment(eps, xb, yb, xc, yc);
      }
    if ((za == Zero) + (zb == Zero) + (zc == Zero) >= 2) return;

    if (za >= Zero && zb >= Zero && zc >= Zero) return;
    if (za <= Zero && zb <= Zero && zc <= Zero) return;

    if ((zc == Zero) || ((zb > Zero) != (zc > Zero)))
      {
        double t; Float s;
        t = xa; xa = xc; xc = t;
        t = ya; ya = yc; yc = t;
        s = za; za = zc; zc = s;
      }
    /* Now zc != 0. */
    /* Also za == 0, or zb == 0, or (sign(za) != sign(zb)) */
    if ((zb == Zero) || ((zb > Zero) != (zc > Zero)))
      {
        double t; Float s;
        t = xa; xa = xb; xb = t;
        t = ya; ya = yb; yb = t;
        s = za; za = zb; zb = s;
      }
    /* Now zb != 0, zc != 0. */
    /* Also za = 0, or (sign(za) != sign(zb) and sign(za) != sign(zc)) */
    if (za != Zero)
      { 
        double sab = zb/(zb - za);
	double sba = za/(za - zb);
	double xab = xa*sab + xb*sba;
	double yab = ya*sab + yb*sba;
	double sac = zc/(zc - za);
	double sca = za/(za - zc);
	double xac = xa*sac + xc*sca;
	double yac = ya*sac + yc*sca;
	epswr_segment(eps, xab, yab, xac, yac);
      }
    else 
      {
	double sbc = zc/(zc - zb);
	double scb = zb/(zb - zc);
	double xbc = xb*sbc + xc*scb;
	double ybc = yb*sbc + yc*scb;
        affirm((zb > 0) != (zc > 0), "fn2_zf_flt_in_triangle: sign confusion");
	epswr_segment(eps, xa, ya, xbc, ybc);
      }
  }

