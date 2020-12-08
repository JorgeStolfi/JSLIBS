/* See pswr_iso.h */
/* Last edited on 2009-08-30 16:11:12 by stolfi */

#define _GNU_SOURCE
#include <math.h>
/* #include <stdio.h> */
/* #include <stdlib.h> */
/* #include <limits.h> */
#include <fpu_control.h>

#include <affirm.h>

#include <pswr_iso.h>

/* ISOLINES */

#define MaxLevels (100000)
  /* Limit to avoid excessive level table allocation. */

#define MaxIsolinesInTriangle (20)
#define MaxIsolinesInRange (20)
  /* Limits to avoid excessive plotting for bad choices of {vStep}. */

static bool_t pswr_isolines_complained = FALSE;
  /* To supress repeated complaints about too many isolines. */

/* LEVELS */

#define DOUBLE_PREC \
  fpu_control_t cw, dw; \
  _FPU_GETCW(cw); \
  dw = (_FPU_IEEE & ~_FPU_EXTENDED & ~_FPU_SINGLE) | _FPU_DOUBLE; \
  _FPU_SETCW(dw)
  
#define RESTORE_FPU \
  _FPU_SETCW(cw)

double pswr_level(double vStart, double vStep, int k)
  { /* Will those #@!$ hackers keep their pizza-stained hands off this code? */
    /* I said DOUBLE, not EXTENDED, dammit!!!! */ 
    DOUBLE_PREC;
    double z = vStart + k*vStep;
    RESTORE_FPU;
    return z;
  }

int pswr_inf_isoline(double vStart, double vStep, double z)
  { 
    affirm(vStep > 0.0, "bad vStep"); 
    DOUBLE_PREC;
    int k = (int)floor((z - vStart)/vStep);
    /* Correct possible roundoff errors: */
    double fk = pswr_level(vStart, vStep, k), fm;
    while ((z < fk) && (fk > (fm = pswr_level(vStart, vStep, k-1)))) 
      { k--; fk = fm; }
    while ((z >= (fm = pswr_level(vStart, vStep, k+1))) && (fm > fk)) 
      { k++; fk = fm; }
    RESTORE_FPU;
    return k;
  }
  
int pswr_sup_isoline(double vStart, double vStep, double z)
  { 
    affirm(vStep > 0.0, "bad vStep"); 
    DOUBLE_PREC;
    int k = (int)ceil((z - vStart)/vStep);
    /* Correct possible roundoff errors: */
    double fk = pswr_level(vStart, vStep, k), fm;
    while ((z > fk) && (fk < (fm = pswr_level(vStart, vStep, k+1)))) 
      { k++; fk = fm; }
    while ((z <= (fm = pswr_level(vStart, vStep, k-1))) && (fm < fk)) 
      { k--; fk = fm; }
    RESTORE_FPU;
    return k;
  }

/* DRAWING ISOLINES */

void pswr_isolines_in_triangle
  ( PSStream *ps,
    double xa, double ya, double fa,
    double xb, double yb, double fb,
    double xc, double yc, double fc,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax        /* Maximum isoline index. */
  )
  { if (kMin > kMax) { return; }

    /* Find the function range {[fMin _ fMax]} in triangle {a,b,c}. */
    double fMin = (fa < fb ? fa : fb);
    fMin = (fMin < fc ? fMin : fc);
    double fMax = (fa > fb ? fa : fb);
    fMax = (fMax > fc ? fMax : fc);
    
    /* Find indices of isolines that intercept the triangle: */
    /* Round values slightly outwards so that we don't lose isolines. */
    double eps = 1.0e-13 * fmin(vStep, fabs(fMax - fMin));
    int iMin = pswr_sup_isoline(vStart, vStep, fMin - eps);
    if (iMin > kMax) { return; }
    if (iMin < kMin) { iMin = kMin; }
    
    int iMax = pswr_inf_isoline(vStart, vStep, fMax + eps);
    if (iMax < kMin) { return; }
    if (iMax > kMax) { iMax = kMax; }

    if (iMin > iMax) { return; }
         
    /* Sanity check to avoid excessive plotting: */
    int nIsolines = iMax - iMin + 1;
    if (nIsolines > MaxIsolinesInTriangle)
      { if (! pswr_isolines_complained) 
          { fprintf(stderr, "** too many isolines in triangle");
            fprintf(stderr, " fMin = %g  fMax = %g", fMin, fMax);
            pswr_isolines_complained = TRUE;
          }
        int excess = MaxIsolinesInTriangle - nIsolines;
        iMin = iMin + excess/2;
        iMax = iMin + MaxIsolinesInTriangle - 1;
        affirm((iMin >= kMin) && (iMax <= kMax), "isoline bug");
      }

    /* fprintf(stderr, "vStart = %.15e vStep = %.15e kMin = %d kMax = %d\n", vStart, vStep, kMin, kMax); */
    /* fprintf(stderr, "fMin = %.15e fMax = %.15e iMin = %d iMax = %d\n", fMin, fMax, iMin, iMax);  */

    /* Permute {a,b,c} so that {fa <= fb <= fc}: */
    pswr_sort_triangle(&xa, &ya, &fa,  &xb, &yb, &fb,  &xc, &yc, &fc);

    /* Plot the isolines: */
    int k;
    for (k = iMin; k <= iMax; k++)
      { double vk = pswr_level(vStart, vStep, k);
        if ((fa == vk) && (fc == vk))
          { /* Degenerate triangle -- all at the same level {vk}. */
            pswr_segment(ps, xa, ya, xb, yb);
            pswr_segment(ps, xb, yb, xc, yc);
            pswr_segment(ps, xc, yc, xa, ya);
          }
        else if ((fa <= vk) && (vk <= fc))
          { double xu, yu, xv, yv;
            pswr_compute_zero_line_in_triangle
              ( xa, ya, fa - vk,
                xb, yb, fb - vk,
                xc, yc, fc - vk,
                &xu, &yu, &xv, &yv
              );
            if ((xu != xv) || (yu != yv))
              { pswr_segment(ps, xu, yu, xv, yv);}
          }
        else 
          { /* Isoline outside {f} range: probable rounding problem, skip it. */ }
      }
  }

void pswr_isolines_in_quadrilateral
  ( PSStream *ps,
    double x00, double y00, double f00,
    double x01, double y01, double f01,
    double x10, double y10, double f10,
    double x11, double y11, double f11,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax        /* Maximum isoline index. */
  )
  { 
    if (kMin > kMax) { return; }
    
    /* Check function range in quadrilateral, skip if no isolines: */
    { double fMin = (f00 < f01 ? f00 : f01);
      fMin = (fMin < f10 ? fMin : f10);
      fMin = (fMin < f11 ? fMin : f11);
      if (fMin > pswr_level(vStart, vStep, kMax)) { return; }
    }
    { double fMax = (f00 > f01 ? f00 : f01);
      fMax = (fMax > f10 ? fMax : f10);
      fMax = (fMax > f11 ? fMax : f11);
      if (fMax < pswr_level(vStart, vStep, kMin)) { return; }
    }
    
    /* Compute coords and estimated value at center of quadrilateral: */
    double xb = (x00 + x01 + x10 + x11)/4.0;
    double yb = (y00 + y01 + y10 + y11)/4.0;
    double fb = (f00 + f01 + f10 + f11)/4.0;
    
    /* Paint (side,center) triangles: */
    pswr_isolines_in_triangle
      ( ps, x00, y00, f00, x10, y10, f10, xb, yb, fb, vStart, vStep, kMin, kMax );
    pswr_isolines_in_triangle
      ( ps, x10, y10, f10, x11, y11, f11, xb, yb, fb, vStart, vStep, kMin, kMax );
    pswr_isolines_in_triangle
      ( ps, x11, y11, f11, x01, y01, f01, xb, yb, fb, vStart, vStep, kMin, kMax );
    pswr_isolines_in_triangle
      ( ps, x01, y01, f01, x00, y00, f00, xb, yb, fb, vStart, vStep, kMin, kMax );
  }

/* PAINTING COLOR BANDS */
  
void pswr_compute_band_indices
  ( double vStart, 
    double vStep, 
    double fMin, 
    int kMin, 
    int *iMin, 
    double *zMin, 
    double fMax, 
    int kMax, 
    int *iMax, 
    double *zMax 
  )
  {
    /* This procedure is tricky to write because the hackers who wrote */
    /* GCC were rather sloppy about "double" versus "extended double" format. */
    /* We must hack the FPU control word to use strict double precision... */
    DOUBLE_PREC;
    if (kMin > kMax) 
      { /* There are no isolines, only one color band: */
        *iMin = *iMax = kMin;
        *zMin = *zMax = pswr_level(vStart, vStep, *iMin);
        RESTORE_FPU; return;
      }
    *zMin = vStart + kMin*vStep;  /* Lowest isoline level. */
    *zMax = vStart + kMax*vStep;  /* Highest isoline level. */
    if (fMin >= *zMax) 
      { /* Triangle is entirely above highest isoline: */
        *iMin = *iMax = kMax + 1;
        *zMin = *zMax = pswr_level(vStart, vStep, *iMin);
      }
    else if (fMax <= *zMin) 
      { /* Triangle is entirely below lowest isoline: */
        *iMin = *iMax = kMin;
        *zMin = *zMax = pswr_level(vStart, vStep, *iMin);
      }
    else
      { /* Triangle intersects some of the bands {kMin..kMax}: */
        /* Get highest level {iMax} in {kMin..kMax+1} s.t. {z[iMax-1] <= fMax}: */
        *iMax = pswr_sup_isoline(vStart, vStep, fMax);
        *zMax = pswr_level(vStart, vStep, *iMax);
        affirm(pswr_level(vStart, vStep, (*iMax)-1) < fMax, "iMax bug 1");
        if (*iMax < kMin)
          { *iMax = kMin; *zMax = pswr_level(vStart, vStep, *iMax); }
        else if (*iMax > kMax) 
          { *iMax = kMax + 1; *zMax = +INFINITY; }
        affirm(fMax <= *zMax, "iMax bug 1");
        
        if (fMin == *zMax) 
          { /* Entire triangle coincides with isoline {iMax}: */
            iMin = iMax;
            *zMin = *zMax;
            RESTORE_FPU; return;
          }

        /* Get lowest isoline {iMin} in {kMin..kMax+1} s.t. {fMin < z[iMin]}: */
        *iMin = pswr_inf_isoline(vStart, vStep, fMin) + 1;
        *zMin = pswr_level(vStart, vStep, *iMin);
        affirm(pswr_level(vStart, vStep, (*iMin)-1) <= fMin, "iMin bug 1");
        if (*iMin < kMin) 
          { *iMin = kMin; *zMin = pswr_level(vStart, vStep, *iMin); }
        else if (*iMin > kMax) 
          { *iMin = kMax + 1; *zMin = INFINITY;  }
        affirm(fMin < *zMin, "iMin bug 2");
        
        affirm(*iMin <= *iMax, "iMin,iMax bug");
      }
    RESTORE_FPU;
  }

void pswr_bands_in_triangle
  ( PSStream *ps,
    double xa, double ya, double fa,
    double xb, double yb, double fb,
    double xc, double yc, double fc,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax,       /* Maximum isoline index. */
    double *R, double *G, double *B
  )
  { /* Find function range in triangle: */
    double fMin = (fa < fb ? fa : fb);
    fMin = (fMin < fc ? fMin : fc);
    double fMax = (fa > fb ? fa : fb);
    fMax = (fMax > fc ? fMax : fc);
    
    /* Find range {iMin .. iMax} of color band indices: */
    int iMin, iMax;
    double zMin, zMax;
    pswr_compute_band_indices
      ( vStart, vStep, 
        fMin, kMin, &iMin, &zMin, 
        fMax, kMax, &iMax, &zMax 
      );

    affirm(iMin <= iMax, "bad bands");

    /* Paint bands: */
    if ((iMin == iMax) || (fMin == fMax))
      { /* Entire triangle lies between levels {iMin-1} (incl) and {iMin} (incl): */
        double Rk = R[iMin - kMin], Gk = G[iMin - kMin], Bk = B[iMin - kMin];
        pswr_set_fill_color(ps, Rk,Gk,Bk);
        pswr_triangle(ps, xa, ya, xb, yb, xc, yc, TRUE, FALSE);
      }
    else
      { /* Triangle has nontrivial intersection with bands {iMin..iMax}: */ 
        /* Permute {a,b,c} so that {fa <= fb <= fc}: */
        pswr_sort_triangle(&xa, &ya, &fa,  &xb, &yb, &fb,  &xc, &yc, &fc);
        affirm(fa == fMin, "bad sort");
        affirm(fc == fMax, "bad sort");
        
        /* Break the triangle {a,b,c} into slices and paint them: */
        double zLo = fa;          /* Value at lower isoline. */
        double xr = xa, yr = ya;  /* Lower point on polyg {a--b--c}. */
        double xs = xa, ys = ya;  /* Lower point on seg {a--c}. */
        bool_t sameLo = TRUE;       /* Do {r} and {s} coincide? */
        int k;
        for (k = iMin; k <= iMax; k++)
          { double zHi;     /* Value at higher isoline. */
            double xu, yu;  /* Higher point on polyg {a--b--c}. */
            double xv, yv;  /* Higher point on seg {a--c}. */
            bool_t sameHi;    /* Do {u} and {v} coincide? */
            
            /* The slice between isolines {k-1} and {k} is defined by:
              two lower points {r=(xr,yr)} on {a--b--c} and
              {s=(xs,ys)} on {a--c}; two higher points {u=(xu,yu)} on
              {a--b--c} and {v=(xv,yv)} on {a--c}; and possibly the
              point {b}, if {zLo < fb < zHi}. The convex region
              delimited by those points lies between isolines {k-1}
              and {k}, where the function values are {zLo} and {zHi},
              respectively. */

            /* Compute the slice's upper edge {u--v}: */
            zHi = (k >= iMax ? zMax : pswr_level(vStart, vStep, k));
            if (zHi >= fc)
              { /* Last band (or whole triangle at this level): */
                xv = xu = xc; yv = yu = yc;
                sameHi = TRUE;
              }
            else if (zHi <= fa)
              { /* Band is empty: */
                xv = xu = xa; yv = yu = ya;
                sameHi = TRUE;
              }
            else 
              { /* fprintf(stderr, "vStart = %.15e vStep = %.15e\n", vStart, vStep); */
                /* fprintf(stderr, "k = %d kMin = %d iMin = %d iMax = %d kMax+1 = %d\n", k, kMin, iMin, iMax, kMax+1); */
                /* fprintf(stderr, "zHi = %.15e fa = %.15e fc = %.15e\n", zHi, fa, fc); */
                pswr_compute_zero_line_in_triangle
                  ( xa, ya, fa - zHi,
                    xb, yb, fb - zHi,
                    xc, yc, fc - zHi,
                    &xu, &yu, &xv, &yv
                  );
                sameHi = FALSE;
              }
            /* Paint the slice: */
            double Rk = R[k - kMin], Gk = G[k - kMin], Bk = B[k - kMin];
            if (Rk >= 0.0) 
              { /* Split the slice into trianglets, and paint them: */
                pswr_set_fill_color(ps, Rk,Gk,Bk);

                if ((zLo > fb) || (zHi < fb))
                  { /* Region is the trapezoid {s,r,u,v}. */
                    if (sameLo && sameHi)
                      { /* Nothing to do */ }
                    else if (sameLo)
                      { pswr_triangle(ps, xs, ys, xu, yu, xv, yv, TRUE, FALSE); }
                    else if (sameHi)
                      { pswr_triangle(ps, xs, ys, xr, yr, xu, yu, TRUE, FALSE); }
                    else
                      { pswr_triangle(ps, xs, ys, xu, yu, xv, yv, TRUE, FALSE);
                        pswr_triangle(ps, xu, yu, xs, ys, xr, yr, TRUE, FALSE);
                      }
                  }
                else
                  { /* Region is the convex pentagon {s,r,b,u,v}. */
                    if (sameLo && sameHi)
                      { pswr_triangle(ps, xs, ys, xb, yb, xv, yv, TRUE, FALSE); }
                    else if (sameLo)
                      { pswr_triangle(ps, xs, ys, xb, yb, xu, yu, TRUE, FALSE);
                        pswr_triangle(ps, xs, ys, xu, yu, xv, yv, TRUE, FALSE);
                      }
                    else if (sameHi)
                      { pswr_triangle(ps, xs, ys, xr, yr, xu, yu, TRUE, FALSE);
                        pswr_triangle(ps, xr, yr, xb, yb, xu, yu, TRUE, FALSE);
                      }
                    else
                      { pswr_triangle(ps, xr, yr, xb, yb, xs, ys, TRUE, FALSE);
                        pswr_triangle(ps, xs, ys, xb, yb, xv, yv, TRUE, FALSE);
                        pswr_triangle(ps, xv, yv, xb, yb, xu, yu, TRUE, FALSE);
                      }
                  }
              }
            /* Prepare for next slice: */
            xr = xu; yr = yu;
            xs = xv; ys = yv;
            zLo = zHi; 
            sameLo = sameHi;
          }
      }
  }

void pswr_bands_in_quadrilateral
  ( PSStream *ps,        /* Plot file. */
    double x00, double y00, double f00,
    double x01, double y01, double f01,
    double x10, double y10, double f10,
    double x11, double y11, double f11,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax,       /* Maximum isoline index. */
    double *R, double *G, double *B
  )
  { /* Compute coords and estimated value at center of quadrilateral: */
    double xb = (x00 + x01 + x10 + x11)/4.0;
    double yb = (y00 + y01 + y10 + y11)/4.0;
    double fb = (f00 + f01 + f10 + f11)/4.0;
    
    /* Paint (side,center) triangles: */
    pswr_bands_in_triangle
      ( ps, x00,y00,f00, x10,y10,f10, xb,yb,fb, vStart, vStep, kMin, kMax, R,G,B );
    pswr_bands_in_triangle
      ( ps, x10,y10,f10, x11,y11,f11, xb,yb,fb, vStart, vStep, kMin, kMax, R,G,B );
    pswr_bands_in_triangle
      ( ps, x11,y11,f11, x01,y01,f01, xb,yb,fb, vStart, vStep, kMin, kMax, R,G,B );
    pswr_bands_in_triangle
      ( ps, x01,y01,f01, x00,y00,f00, xb,yb,fb, vStart, vStep, kMin, kMax, R,G,B );
  }       

void pswr_compute_zero_line_in_triangle
  ( double xa, double ya, double fa,
    double xb, double yb, double fb,
    double xc, double yc, double fc,
    double *xu, double *yu,
    double *xv, double *yv
  )
  { affirm(fa <= fb , "sort error(fa,*fb)");
    affirm(fb <= fc , "sort error(fb,*fc)");
    affirm(fa < fc , "degenerate zero line");
    
    /* If the zero set is empty or a single point, skip: */
    if ((fa >= 0.0) && (fb > 0.0)) 
      { *xu = *xv = xa; *yu = *yv = ya; }
    else if ((fc <= 0.0) && (fb < 0.0)) 
      { *xu = *xv = xc; *yu = *yv = yc; }
    else
      { /* Find {f = 0} point {u} on polygonal {a--b--c}: */
        if (fb < 0.0) 
          { double s_bc = fc/(fc - fb);
            double s_cb = fb/(fb - fc);
            *xu = xb*s_bc + xc*s_cb;
            *yu = yb*s_bc + yc*s_cb;
          }
        else if (fb > 0.0) 
          { double s_ba = fa/(fa - fb);
            double s_ab = fb/(fb - fa);
            *xu = xb*s_ba + xa*s_ab;
            *yu = yb*s_ba + ya*s_ab;
          }
        else
          { *xu = xb; *yu = yb; }

        /* Find {f = 0} point {v} on side {a--c}: */
        if (fa == 0.0)
          { *xv = xa; *yv = ya; }
        else if (fc == 0.0)
          { *xv = xc; *yv = yc; }
        else
          { double s_ca = fa/(fa - fc);
            double s_ac = fc/(fc - fa);
            *xv = xc*s_ca + xa*s_ac;
            *yv = yc*s_ca + ya*s_ac;
          }
      }
  }

void pswr_sort_triangle
  ( double *xa, double *ya, double *fa,
    double *xb, double *yb, double *fb,
    double *xc, double *yc, double *fc
  )
  /* Permute corners so that {fc <= fd <= fb}: */
  {
    if (*fa > *fb)
      { double temp;
        temp = *xb; *xb = *xa; *xa = temp;
        temp = *yb; *yb = *ya; *ya = temp;
        temp = *fb; *fb = *fa; *fa = temp;
      }
    if (*fb > *fc)
      { double temp;
        temp= *xb; *xb = *xc; *xc = temp;
        temp= *yb; *yb = *yc; *yc = temp;
        temp= *fb; *fb = *fc; *fc = temp;
      }
    if (*fa > *fb)
      { double temp;
        temp = *xb; *xb = *xa; *xa = temp;
        temp = *yb; *yb = *ya; *ya = temp;
        temp = *fb; *fb = *fa; *fa = temp;
      }
  }    


