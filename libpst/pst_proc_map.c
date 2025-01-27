/* See pst_proc_map.h */
/* Last edited on 2025-01-26 19:42:28 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

#include <bool.h>
#include <vec.h>
#include <float_image.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <r2.h>
#include <r2x2.h>
#include <r3.h>

#include <pst_proc_map.h>
#include <pst_proc_map_fractal.h>
#include <pst_normal_map.h>

/* INTERNAL PROTOTYPES */

void pst_proc_map_compute_average_gradient_and_weight
  ( double sum_wx_dzx,
    double sum_wx,
    double sum_wy_dzy,
    double sum_wy,
    bool_t ok,
    double maxGrad,
    bool_t debug,
    r2_t *dz_P,
    double *wdz_P
  );
  /* Computes the average gradient {dz} and its reliability weight
    {wdz}, given, for each axis, the sums of weights {sum_wx,xum_wy} and
    sums weights times derivatives {sum_wx_dzx,sum_wy_dzy} included in
    the average, and the flag {ok} which should be false if any would-be
    term was invalid or exceeded {maxGrad}. Returns the results in
    {*dz_P} and {*wdz_P}.
    
    The gradient {dz} will be {(sum_wx_dzx/sum_wx, sum_wy_dzysum_wy)}
    and the weight {wdz} will be 1.0 if {ok} is true and both quotients
    are finite, and its norm does not exceed {maxGrad}. Otherwise {dz}
    will be {(NAN,NAN)} and {wdz} will be zero. */

/* IMPLEMENTATIONS */

void pst_proc_map_compute_height
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    uint32_t NS,
    double ws[],
    double xyScale,
    bool_t debug,
    double *z_P,
    double *wz_P
  )
  {
    demand(isfinite(xyScale) && (xyScale > 0), "invalid {xyScale}");
    demand(((NS%2) == 1) && (NS >= 3), "sampoint count {NS} must be odd and at least 3");

    int32_t HS = (int32_t)NS/2;
    double step = 1.0/(HS+1.0)/xyScale; /* Displacement between sampoints. */

    /* Accumulators for average height: */
    double sum_wz = 0;    /* Sum of weight times height for valid sampoints. */
    double sum_w = 0;     /* Sum of weights for valid sampoints. */
    double tot_w = 0;     /* Sum of weights for all sampoints. */

    for (int32_t ky = -HS; ky <= +HS; ky++)
      { for (int32_t kx = -HS; kx <= +HS; kx++)
          { /* Generate a sampoint {(xk,yk)} around {p}: */
            double xk = p.c[0] + ((double)kx)*step;
            double yk = p.c[1] + ((double)ky)*step;
            r2_t pk = (r2_t){{ xk, yk }};
            if (debug) { fprintf(stderr, "  %3d %3d pk = ( %12.8f %12.8f )", kx, ky, xk, yk); }

            /* Get height value: */
            double zk;
            func(pk, &zk, NULL);

            /* Get the sampoint weight {wk}: */
            double wk = ws[HS+kx]*ws[HS+ky];
            tot_w += wk;

            if (isfinite(zk))
              { /* Accumulate height value: */
                if (debug) { fprintf(stderr, "  zk = %24.16e  wk = %24.16e\n", zk, wk); }
                sum_wz += wk*zk;
                sum_w += wk;
              }
         }
     }

    /* Compute mean height: */
    double z_avg = (sum_w == 0 ? NAN: sum_wz/sum_w);
    double wz_avg = (isnan(z_avg) ? 0.0 : (sum_w/tot_w));
    if (debug) { fprintf(stderr, "  avg z = %24.16e wz = %24.16e\n", z_avg, wz_avg); }

    /* Return results: */
    if (z_P != NULL) { (*z_P) = z_avg; }
    if (wz_P != NULL) { (*wz_P) = wz_avg; }
  }

void pst_proc_map_compute_numeric_gradient
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    uint32_t NS,
    double ws[],
    double xyScale,
    double maxGrad,
    bool_t debug,
    r2_t *dnz_P,
    double *wdnz_P
  )
  {
    demand(isfinite(xyScale) && (xyScale > 0), "invalid {xyScale}");
    demand((! isnan(maxGrad)) && (maxGrad > 0), "invalid maxGrad");
    demand(((NS%2) == 1) && (NS >= 3), "sampoint count {NS} must be odd and at least 3");

    int32_t HS = (int32_t)NS/2;
    double step = 1.0/(HS+1.0)/xyScale; /* Displacement between sampoiints. */

    /* Accumulators for sampoints that are not {NAN}: */
    r2_t sum_wn_dz = (r2_t){{ 0, 0 }};
    r2_t sum_wn = (r2_t){{ 0, 0 }};
    
    bool_t ok = TRUE; /* Set to false if any sampoint or pair is invalid. */

    double zprev[NS]; /* Height values on previous sampoint row. */

    for (int32_t ky = -HS; ok && (ky <= +HS); ky++)
      { for (int32_t kx = -HS; ok && (kx <= +HS); kx++)
          { /* Generate a sampoint {(xk,yk)} around {p}: */
            double xk = p.c[0] + ((double)kx)*step;
            double yk = p.c[1] + ((double)ky)*step;
            r2_t pk = (r2_t){{ xk, yk }};
            if (debug) { fprintf(stderr, "  %3d %3d pk = ( %12.8f %12.8f )\n", kx, ky, xk, yk); }

            /* Get height value and analytic derivatives: */
            double zk;
            func(pk, &zk, NULL);

            /* Get the weights of the pairs: */
            double wnxk = (HS+kx > 0 ? 0.5*(ws[HS+kx] + ws[HS+kx-1])*ws[HS+ky] : 0.0);
            double wnyk = (HS+ky > 0 ? 0.5*(ws[HS+ky] + ws[HS+ky-1])*ws[HS+kx] : 0.0);

            /* Accumulate numeric gradient in {sum_wn_dz,sum_wn}: */
            if (! isfinite(zk))
              { ok = FALSE; }
            else
              { if (HS+kx > 0)
                  { /* Compute and accumulate numeric {X} derivative: */
                    double zpxk = zprev[HS+kx-1];
                    assert(isfinite(zpxk)); /* Should have stopped if not. */
                    double dnzxk = (zk - zpxk)/step;
                    if (debug) { fprintf(stderr, "  dnzxk = %24.16e  wnxk = %24.16e\n", dnzxk, wnxk); }
                    if (fabs(dnzxk) > maxGrad) 
                      { ok = FALSE; }
                    else
                      { sum_wn_dz.c[0] += wnxk*dnzxk;
                        sum_wn.c[0] += wnxk;
                      }
                  }
                if (HS+ky > 0)
                  { /* Compute and accumulate numeric {Y} derivative: */
                    double zpyk = zprev[HS+kx];
                    assert(isfinite(zpyk)); /* Should have stopped if not. */
                    double dnzyk = (zk - zpyk)/step;
                    if (debug) { fprintf(stderr, "  dnzyk = %24.16e  wnyk = %24.16e\n", dnzyk, wnyk); }
                    if (fabs(dnzyk) > maxGrad)
                      { ok = FALSE; }
                    else
                      { sum_wn_dz.c[1] += wnyk*dnzyk;
                        sum_wn.c[1] += wnyk;
                      }
                  }
              }

            if (debug) { fprintf(stderr, "\n"); }
            zprev[HS+kx] = zk;
          }
      }
    r2_t dnz_avg; double wdnz_avg;
    pst_proc_map_compute_average_gradient_and_weight
      ( sum_wn_dz.c[0], sum_wn.c[0],
        sum_wn_dz.c[1], sum_wn.c[1],
        ok, maxGrad, debug,
        &dnz_avg, &wdnz_avg
      );
    if (debug) { fprintf(stderr, "  avg dnzx = %24.16e dnzy = %24.16e\n", dnz_avg.c[0], dnz_avg.c[1]); }

    /* Return results: */
    if ( dnz_P != NULL) { (*dnz_P) = dnz_avg; }
    if ( wdnz_P != NULL) { (*wdnz_P) = wdnz_avg; }
  }

void pst_proc_map_compute_analytic_gradient
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    uint32_t NS,
    double ws[],
    double xyScale,
    double maxGrad,
    bool_t debug,
    r2_t *daz_P,
    double *wdaz_P
  )
  {
    demand(isfinite(xyScale) && (xyScale > 0), "invalid {xyScale}");
    demand((! isnan(maxGrad)) && (maxGrad > 0), "invalid maxGrad");
    demand(((NS%2) == 1) && (NS >= 3), "sampoint count {NS} must be odd and at least 3");

    int32_t HS = (int32_t)NS/2;
    double step = 1.0/(HS+1.0)/xyScale; /* Displacement between sampoiints. */

    /* Accumulators for average analytic gradient: */
    r2_t sum_wa_dz = (r2_t){{ 0, 0 }}; /* Sum of weight × derivative for valid sampoints. */
    double sum_wa = 0; /* Sum of weights for valid sampoints. */

    bool_t ok = TRUE; /* Set to false if any sampoint or pair is invalid. */

    for (int32_t ky = -HS; ok && (ky <= +HS); ky++)
      { for (int32_t kx = -HS; ok && (kx <= +HS); kx++)
          { /* Generate a sampoint {(xk,yk)} around {p}: */
            double xk = p.c[0] + ((double)kx)*step;
            double yk = p.c[1] + ((double)ky)*step;
            r2_t pk = (r2_t){{ xk, yk }};
            if (debug) { fprintf(stderr, "  %3d %3d pk = ( %12.8f %12.8f )", kx, ky, xk, yk); }

            /* Get height value and analytic derivatives: */
            r2_t dazk;
            func(pk, NULL, &dazk);

            /* Get the sampoint weight {wk}: */
            double wak = ws[HS+kx]*ws[HS+ky];

            if (isfinite(dazk.c[0]) && isfinite(dazk.c[1]))
              { if (debug) { fprintf(stderr, "  dazk = ( %24.16e %24.16e ) wak = %24.16e", dazk.c[0], dazk.c[1], wak); }
                /* Check max gradient condition: */
                double dazmk = hypot(dazk.c[0], dazk.c[1]);
                if (dazmk > maxGrad)
                  { ok = FALSE; }
                else
                  { /* Accumulate analytic slopes in {sum_wa_dz,sum_wa}: */
                    sum_wa_dz.c[0] += wak*dazk.c[0];
                    sum_wa_dz.c[1] += wak*dazk.c[1];
                    sum_wa += wak;
                  }
              }
            else
              { ok = FALSE; }

            if (debug) { fprintf(stderr, "\n"); }
         }
     }

    r2_t daz_avg; double wdaz_avg;
    pst_proc_map_compute_average_gradient_and_weight
      ( sum_wa_dz.c[0], sum_wa,
        sum_wa_dz.c[1], sum_wa,
        ok, maxGrad, debug,
        &daz_avg, &wdaz_avg
      );
    if (debug) { fprintf(stderr, "  avg dazx = %24.16e dazy = %24.16e\n", daz_avg.c[0], daz_avg.c[1]); }

    /* Return results: */
    if ( daz_P != NULL) { (*daz_P) = daz_avg; }
    if ( wdaz_P != NULL) { (*wdaz_P) = wdaz_avg; }
  }

void pst_proc_map_compute_average_gradient_and_weight
  ( double sum_wx_dzx,
    double sum_wx,
    double sum_wy_dzy,
    double sum_wy,
    bool_t ok,
    double maxGrad,
    bool_t debug,
    r2_t *dz_P,
    double *wdz_P
  )
  { /* Compute mean numeric gradient: */
    double dzx_avg = (sum_wx == 0 ? NAN : sum_wx_dzx/sum_wx);
    double dzy_avg = (sum_wy == 0 ? NAN : sum_wy_dzy/sum_wy);
    r2_t dz_avg;
    if ((! ok) || (! isfinite(dzx_avg)) || (! isfinite(dzy_avg)))
      { dz_avg = (r2_t){{ NAN, NAN }}; }
    else
      { /* Check again the {maxGrad} condition, just in case: */
        double dzm = hypot(dzx_avg, dzy_avg);
        if (dzm > maxGrad)
          { dz_avg = (r2_t){{ NAN, NAN }}; }
        else
          { dz_avg = (r2_t){{ dzx_avg, dzy_avg }}; }
      }

    /* Compute overal reliability weight: */
    double wdz_avg;
    if ((! ok) || isnan(dz_avg.c[0]) || isnan(dz_avg.c[1]))
      { wdz_avg = 0.0; }
    else
      { wdz_avg = 1.0; }

    /* Return results: */
    (*dz_P) = dz_avg;
    (*wdz_P) = wdz_avg;
  }

float_image_t* pst_proc_map_make_height_map
  ( pst_proc_map_zfunc_t *func,
    int32_t NX,
    int32_t NY,
    uint32_t NS,
    double ws[]
  )
  {
    int32_t xDebug = -1, yDebug = -1;

    int32_t NC = 2;
    float_image_t *IZ = float_image_new(NC, NX+1, NY+1);

    int32_t NMIN = (NX < NY ? NX : NY);  /* Smallest dimension of slope map. */
    double xyScale = ((double)NMIN)/2.0;  /* Number of grid pixels per {func} domain unit. */
    /* Center of image domain (in pixels, from bottom left corner): */
    double cX = 0.5*NX;
    double cY = 0.5*NY;

    /* Compute height at each grid CORNER: */
    for (int32_t y = 0; y <= NY; y++)
      { for (int32_t x = 0; x <= NX; x++)
          { bool_t debug = ((x == xDebug) && (y == yDebug));

            /* Compute function domain coordinates {(xp,yp)} of grid corner: */
            double xp = (x - cX)/xyScale;
            double yp = (y - cY)/xyScale;
            r2_t p = (r2_t){{ xp, yp }};

            /* Compute height at grid corner: */
            double z, wz;
            if (debug) { fprintf(stderr, "=== debugging {pst_proc_map_compute_height} for pixel [%d,%d]\n", x, y); }
            pst_proc_map_compute_height(func, p, NS, ws, xyScale, debug, &z, &wz);
            if (debug) { fprintf(stderr, "=== end debugging\n"); }

            /* Scale height from {func} units to pixels: */
            float_image_set_sample(IZ, 0, x, y, (float)(z*xyScale));
            float_image_set_sample(IZ, 1, x, y, (float)wz);
          }
      }
    return IZ;
  }

float_image_t* pst_proc_map_make_slope_map
  ( pst_proc_map_zfunc_t *func,
    int32_t NX,
    int32_t NY,
    uint32_t NS,
    double ws[],
    bool_t numGrad,
    double maxGrad,
    double maxGDiff
  )
  {
    int32_t xDebug = 197, yDebug = 96;

    int32_t NC = 3;
    float_image_t *IG = float_image_new(NC, NX, NY);

    int32_t NMIN = (NX < NY ? NX : NY);   /* Smallest dimension of image. */
    double xyScale = ((double)NMIN)/2.0;  /* Number of grid pixels per {func} domain unit. */

    /* Center of image domain (in pixels, from bottom left corner): */
    double cX = 0.5*NX;
    double cY = 0.5*NY;

    /* Compute gradient at center of each grid pixel: */
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { bool_t debug = ((x == xDebug) && (y == yDebug));

            /* Compute function domain coordinates {(xp,yp)} of pixel center: */
            double xp = (x + 0.5 - cX)/xyScale;
            double yp = (y + 0.5 - cY)/xyScale;
            r2_t p = (r2_t){{ xp, yp }};

            /* Compute gradient at pixel center: */
            if (debug) { fprintf(stderr, "=== debugging {pst_proc_map_compute_gradient} for pixel [%d,%d]\n", x, y); }
            r2_t dnz, daz; double wdnz, wdaz;
            pst_proc_map_compute_numeric_gradient
              ( func, p, NS, ws, xyScale, maxGrad, debug, &dnz, &wdnz );
            pst_proc_map_compute_analytic_gradient
              ( func, p, NS, ws, xyScale, maxGrad, debug, &daz, &wdaz );

            /* If gradients are discrepant, the analytic one is unreliable: */
            if ((wdaz == 0) || (wdnz == 0))
              { wdaz = 0.0; }
            else
              { double ddz = r2_dist(&dnz, &daz);
                if (ddz > maxGDiff) { wdaz = 0.0; }
              }
            if (wdaz == 0) { daz = (r2_t){{ NAN, NAN }}; }

            if (debug) { fprintf(stderr, "=== end debugging\n"); }
            r2_t *dz = (numGrad ? &dnz : &daz);
            double wdz = (numGrad ? wdnz : wdaz);
            float_image_set_sample(IG, 0, x, y, (float)dz->c[0]);
            float_image_set_sample(IG, 1, x, y, (float)dz->c[1]);
            float_image_set_sample(IG, 2, x, y, (float)wdz);
          }
      }
    return IG;
  }

pst_proc_map_zfunc_props_t pst_proc_map_function_generic(int32_t n)
  {
    pst_proc_map_zfunc_props_t fp;

    auto void dfp(pst_proc_map_zfunc_t *func, char *name, double maxGrad, double maxGDiff);

    switch(n)
      { case  0: dfp(&pst_proc_map_function_00, "zeflat", +INF, +INF); break; /* Constant function (zero gradient). */
        case  1: dfp(&pst_proc_map_function_01, "ramp10", +INF, +INF); break; /* Linear ramp in the X direction. */
        case  2: dfp(&pst_proc_map_function_02, "ramp01", +INF, +INF); break; /* Linear ramp in the Y direction. */
        case  3: dfp(&pst_proc_map_function_03, "ramp11", +INF, +INF); break; /* Linear ramp in the XY diagonal direction. */
        case  4: dfp(&pst_proc_map_function_04, "parabo", +INF, +INF); break; /* Slanted elliptical parabolic hump. */
        case  5: dfp(&pst_proc_map_function_05, "spdom1", +INF, +INF); break; /* Buried sphere. */
        case  6: dfp(&pst_proc_map_function_06, "pyram5", +INF, +INF); break; /* Pentagonal pyramid. */
        case  7: dfp(&pst_proc_map_function_07, "rtcone", +INF, +INF); break; /* Conical mound. */
        case  8: dfp(&pst_proc_map_function_08, "ripple", +INF, +INF); break; /* Circular waves. */
        case  9: dfp(&pst_proc_map_function_09, "sbabel", +INF, +INF); break; /* Tower of Babel with smooth shoulders. */
        case 10: dfp(&pst_proc_map_function_10, "hcliff", 5.00, 0.05); break; /* Diverging parabolic ramps with cliff. */
        case 11: dfp(&pst_proc_map_function_11, "sinw01", +INF, +INF); break; /* Sinusoidal wave, low freq. */
        case 12: dfp(&pst_proc_map_function_12, "cosw01", +INF, +INF); break; /* Co-sinusoidal wave, low freq. */
        case 13: dfp(&pst_proc_map_function_13, "sinw02", +INF, +INF); break; /* Sinusoidal wave, med freq. */
        case 14: dfp(&pst_proc_map_function_14, "cosw02", +INF, +INF); break; /* Co-sinusoidal wave, med freq. */
        case 15: dfp(&pst_proc_map_function_15, "sinw03", +INF, +INF); break; /* Sinusoidal wave, high freq. */
        case 16: dfp(&pst_proc_map_function_16, "cosw03", +INF, +INF); break; /* Co-sinusoidal wave, high freq. */
        case 17: dfp(&pst_proc_map_function_17, "cbabel", 5.00, 0.05); break; /* Tower of Babel with cliff. */
        case 18: dfp(&pst_proc_map_function_18, "cbramp", 5.00, 0.05); break; /* Cubic ramp with cliff on three sides. */
        case 19: dfp(&pst_proc_map_function_19, "holes1", +INF, 0.05); break; /* Buried sphere with many holes in slope map. */
        case 20: dfp(&pst_proc_map_function_20, "cpiece", 5.00, 0.05); break; /* Three nested stages connected by ramps. */
        case 21: dfp(&pst_proc_map_function_21, "cplat3", +INF, +INF); break; /* Three flat round platforms with smooth edges. */
        case 22: dfp(&pst_proc_map_function_22, "plat3a", +INF, +INF); break; /* Triangular wall with smooth edges. */
        case 23: dfp(&pst_proc_map_function_23, "plat5a", +INF, +INF); break; /* Pentagonal platform with smooth edges. */
        case 24: dfp(&pst_proc_map_function_24, "cplatv", 5.00, 0.05); break; /* Circular platform with cliffs. */
        case 25: dfp(&pst_proc_map_function_25, "cplats", +INF, +INF); break; /* Circular platform with smooth edges. */
        case 26: dfp(&pst_proc_map_function_26, "qtbell", +INF, +INF); break; /* Quartic bell. */
        case 27: dfp(&pst_proc_map_function_27, "fracto", +INF, +INF); break; /* Fractalish round montain. */
        default:
          fprintf(stderr, "bad {n} %d\n", n);
          assert(FALSE);
      }
    return fp;

    void dfp(pst_proc_map_zfunc_t *func, char *name, double maxGrad, double maxGDiff)
      { fp.num = n;
        fp.func = func;
        fp.name = name;
        fp.maxGrad = maxGrad;
        fp.maxGDiff = maxGDiff;
      }
  }

void pst_proc_map_function_00(r2_t p, double *z, r2_t *dz)
  { /* "zeflat" - constant function: */
    double zval = 0.50;
    if (z != NULL) { (*z) = zval; }
    if (dz != NULL) { dz->c[0] = dz->c[1] = 0; }
  }

void pst_proc_map_function_01(r2_t p, double *z, r2_t *dz)
  { /* "ramp10" - linear ramp along X: */
    double x = p.c[0], y = p.c[1];
    double Cx = 0.5, Cy = 0.0;
    if (z != NULL) { (*z) = Cx*x + Cy*y; }
    if (dz != NULL) { dz->c[0] = Cx; dz->c[1] = Cy; }
  }

void pst_proc_map_function_02(r2_t p, double *z, r2_t *dz)
  { /* "ramp01" - linear ramp along Y: */
    double x = p.c[0], y = p.c[1];
    double Cx = 0.0, Cy = 0.5;
    if (z != NULL) { (*z) = Cx*x + Cy*y; }
    if (dz != NULL) { dz->c[0] = Cx; dz->c[1] = Cy; }
  }

void pst_proc_map_function_03(r2_t p, double *z, r2_t *dz)
  { /* "ramp11" - linear ramp along diagonal: */
    double x = p.c[0], y = p.c[1];
    double Cx = 1/4.0, Cy = 1/3.0;
    if (z != NULL) { (*z) = Cx*x + Cy*y; }
    if (dz != NULL) { dz->c[0] = Cx; dz->c[1] = Cy; }
  }

void pst_proc_map_function_04(r2_t p, double *z, r2_t *dz)
  { /* "parabo" - slanted parabolic hump: */
    double x = p.c[0], y = p.c[1];
    double A = 0.10;
    double B = 0.05;
    double Cx = 1.0, Cy = 1.5;
    double S = Cx*x + Cy*y;
    double T = -Cy*x + Cx*y;
    if (z != NULL) { (*z) = 1.0 - A*S*S - B*T*T; }
    if (dz != NULL)
      { dz->c[0] = -2.0*(A*S*Cx - B*T*Cy);
        dz->c[1] = -2.0*(A*S*Cy + B*T*Cx);
      }
    }

void pst_proc_map_function_05(r2_t p, double *z, r2_t *dz)
  { /* "spdom1" - buried sphere: */
    double x = p.c[0], y = p.c[1];
    double zMin = 0.05; /* Ground level. */
    /* Initialize {*z,*dz} with the ground plane: */
    if (z != NULL) { (*z) = zMin; }
    if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }
    /* Where the sphere is above the ground, set {*z,*dz} to it: */
    double R = 0.8;     /* Sphere radius. */
    double D2 = x*x + y*y;
    if (D2 < R*R)
      { double S = R*R - D2;
        double dSdx = -2*x;
        double dSdy = -2*y;
        double F = sqrt(S);
        if (F > zMin)
          { double dFdS = 0.5/F;
            if (z != NULL) { (*z) = F; }
            if (dz != NULL)
              { dz->c[0] = dFdS*dSdx;
                dz->c[1] = dFdS*dSdy;
              }
          }
      }
  }

void pst_proc_map_function_06(r2_t p, double *z, r2_t *dz)
  { /* "pyram5" - pentagonal pyramid: */
    double x = p.c[0], y = p.c[1];
    double zMin = 0.1; /* Ground level. */
    /* Initialize {*z,*dz} with the ground plane: */
    if (z != NULL) { (*z) = zMin; }
    if (dz != NULL){ dz->c[0] = dz->c[1] = 0.00; }
    /* Where the pyramid is above the ground, set {*z,*dz} to it: */
    double R = 0.75; /* Pyramid radius at z = 0. */
    double H = 0.80; /* Pyramid height. */
    double D2 = x*x + y*y;
    if (D2 == 0)
      { /* At the tip, assume slope = 0: */
        if (z != NULL) { (*z) = H; }
      }
    else
      { double ang = atan2(y,x);
        double dang = 2*M_PI/5;
        double rang = dang*((int32_t)floor((ang + 2*M_PI)/dang + 0.83) - 0.33);
        double rx = cos(rang), ry = sin(rang);
        double S = rx*x + ry*y;
        if (S < R)
          { double dSdx = rx;
            double dSdy = ry;
            double F = H*(1 - S/R);
            if (F > zMin)
              { double dFdS = -H/R;
                if (z != NULL) { (*z) = F; }
                if (dz != NULL)
                  { dz->c[0] = dFdS*dSdx;
                    dz->c[1] = dFdS*dSdy;
                  }
              }
          }
      }
  }

void pst_proc_map_function_07(r2_t p, double *z, r2_t *dz)
  { /* "rtcone" - conical mound: */
    double x = p.c[0], y = p.c[1];
    double zMin = 0.1; /* Ground level. */
    /* Initialize {*z,*dz} with the ground plane: */
    if (z != NULL) { (*z) = zMin; }
    if (dz != NULL){ dz->c[0] = dz->c[1] = 0.00; }
    double R = 0.75; /* Cone radius at z = 0. */
    double H = 0.80; /* Cone height. */
    double D2 = x*x + y*y;
    if (D2 == 0)
      { /* At the tip, assume slope = 0: */
        if (z != NULL) { (*z) = H; }
      }
    else
      { /* Where the cone is above the ground, set {*z,*dz} to it: */
        if (D2 < R*R)
          { double S = sqrt(D2);
            double dSdx = x/S;
            double dSdy = y/S;
            double F = H*(1 - S/R);
            if (F > zMin)
              { double dFdS = -H/R;
                if (z != NULL) { (*z) = F; }
                if (dz != NULL)
                  { dz->c[0] = dFdS*dSdx;
                    dz->c[1] = dFdS*dSdy;
                  }
              }
          }
      }
  }

void pst_proc_map_function_08(r2_t p, double *z, r2_t *dz)
  { /* "ripple" - circular waves: */
    double x = p.c[0], y = p.c[1];
    double A = 0.450; /* Ripple amplitude (one half of peak-to-peak). */
    double R = 0.075; /* Min ripple radius. */
    double W = 1.5*M_PI/M_LN2;  /* Frequency scale factor. */
    double R2 = R*R;
    double D2 = x*x + y*y + R*R;
    double t = log(D2/R2);
    double dtdx = 2*x/D2;
    double dtdy = 2*y/D2;
    double F = A*cos(W*t);
    double dFdt = -A*W*sin(W*t);
    if (z != NULL) { (*z) = F; }
    if (dz != NULL)
      { dz->c[0] = dFdt*dtdx;
        dz->c[1] = dFdt*dtdy;
      }
  }

void pst_proc_map_function_09(r2_t p, double *z, r2_t *dz)
  { /* "sbabel" - tower of Babel with soulders: */
    double N = 3.0;       /* Number of full turns. */
    double RI = 0.05;     /* Inner radius of top platform (excl. shoulder). */
    double RO = 0.90;     /* Outer radius of tower (excl. shoulder). */
    double EF = 1.0/3.0;  /* Relative width of soft shoulder. */

    pst_proc_map_function_babel(&p, RI, RO, N, EF, z, dz);
  }

void pst_proc_map_function_10(r2_t p, double *z, r2_t *dz)
  { /* "hcliff" - diverging parabolic ramps, with cliff: */
    double x = p.c[0], y = p.c[1];
    double Cx = 1.0, Cy = 1.5;
    double S = Cx*x + Cy*y;
    double T = -Cy*x + Cx*y;
    double A = 0.10;
    if (S < 0)
      { /* Flat region: */
        if (z != NULL) { (*z) = 0; }
        if (dz != NULL)
          { dz->c[0] = 0;
            dz->c[1] = 0;
          }
      }
    else if (T > 0)
      { /* Rising parabolic ramp: */
        if (z != NULL) { (*z) = + A*S*S; }
        if (dz != NULL)
          { dz->c[0] = + 2*A*S*Cx;
            dz->c[1] = + 2*A*S*Cy;
          }
      }
    else
      { /* Descending parabolic ramp: */
        if (z != NULL) { (*z) = - A*S*S; }
        if (dz != NULL)
          { dz->c[0] = - 2*A*S*Cx;
            dz->c[1] = - 2*A*S*Cy;
          }
      }
  }

void pst_proc_map_function_11(r2_t p, double *z, r2_t *dz)
  { /* "sinw01" - low-frequency sinusoid wave: */
    r2_t f = (r2_t){{ 0.25, 0.50 }};
    pst_proc_map_function_wave(&p, &f, 0.0, z, dz);
  }

void pst_proc_map_function_12(r2_t p, double *z, r2_t *dz)
  { /* "cosw01" - low-frequency co-sinusoid wave: */
    r2_t f = (r2_t){{ 0.25, 0.50 }};
    pst_proc_map_function_wave(&p, &f, M_PI/2, z, dz);
  }

void pst_proc_map_function_13(r2_t p, double *z, r2_t *dz)
  { /* "sinw02" - medium-frequency sinusoid wave: */
    r2_t f = (r2_t){{ 3.00, 2.50 }};
    pst_proc_map_function_wave(&p, &f, 0.0, z, dz);
  }

void pst_proc_map_function_14(r2_t p, double *z, r2_t *dz)
  { /* "cosw02" - medium-frequency co-sinusoid wave: */
    r2_t f = (r2_t){{ 3.00, 2.50 }};
    pst_proc_map_function_wave(&p, &f, M_PI/2, z, dz);
  }

void pst_proc_map_function_15(r2_t p, double *z, r2_t *dz)
  { /* "sinw03" - high-frequency sinusoid wave: */
    r2_t f = (r2_t){{ 5.00, 7.00 }};
    pst_proc_map_function_wave(&p, &f, 0.0, z, dz);
  }

void pst_proc_map_function_16(r2_t p, double *z, r2_t *dz)
  { /* "cosw03" - higher-frequency co-sinusoid wave: */
    r2_t f = (r2_t){{ 5.00, 7.00 }};
    pst_proc_map_function_wave(&p, &f, M_PI/2, z, dz);
  }

void pst_proc_map_function_17(r2_t p, double *z, r2_t *dz)
  { /* "cbabel" - tower of Babel with cliff: */
    double N = 2.25;      /* Number of full turns. */
    double RI = 0.05;     /* Inner radius of top platform. */
    double RO = 0.90;     /* Outer radius of tower. */
    double EF = 0.0;      /* Relative width of soft shoulder. */

    pst_proc_map_function_babel(&p, RI, RO, N, EF, z, dz);
  }

void pst_proc_map_function_18(r2_t p, double *z, r2_t *dz)
  { /* "cbramp" - tilted bicubic ramp with cliffs on three sides: */
    double x = p.c[0], y = p.c[1];

    /* Tilt {theta}, half-side {A} and height {H} of ramp: */
    double theta = M_PI/7;
    double Ct = cos(theta);
    double St = sin(theta);
    double A = 0.75/(fabs(Ct) + fabs(St));
    double H = 0.80;

    /* Remap {x,y} to {u,v} tilted by theta: */
    double u = (+ Ct*x + St*y)/A;
    double v = (- St*x + Ct*y)/A;

    double Vmax = 0.75;
    double Umax = 1.50;

    if ((u < -1.00) || (u > +Umax) || (fabs(v) > Vmax))
      { /* Ground floor: */
        if (z != NULL) { (*z) = 0; }
        if (dz != NULL)
          { dz->c[0] = 0;
            dz->c[1] = 0;
          }
      }
    else if (u > +1.00)
      { /* Upper floor: */
        if (z != NULL) { (*z) = H; }
        if (dz != NULL)
          { dz->c[0] = 0;
            dz->c[1] = 0;
          }
      }
    else
      { /* Ramp: */
        if (z != NULL) { (*z) = 0.25*H*(2 + 3*u - u*u*u); }
        if (dz != NULL)
          { double dzdu = 0.75*H*(1 - u*u)/A;
            dz->c[0] = dzdu*Ct;
            dz->c[1] = dzdu*St;
          }
      }
  }

void pst_proc_map_function_19(r2_t p, double *z, r2_t *dz)
  { /* "holes1" - buried sphere with undefined regions and scattered outliers: */

    /* Compute the buried sphere without noise: */
    pst_proc_map_function_05(p, z, dz);

    if (dz == NULL) { return; }

    /* Now choose which places will be bad: */
    bool_t bad = FALSE; /* For now. */
    double x = p.c[0], y = p.c[1];

    /* Make it bad outside the dome: */
    double R = 0.95;    /* Slightly bigger than sphere radius. */
    double D2 = x*x + y*y;
    if (D2 > R*R) { bad = TRUE; }

    /* Make a pizza slice bad: */
    if ((y > 0) && (y < -x/3)) { bad = TRUE; }

    /* Drill some holes in various places: */
    int32_t M = 2;
    double xf = M*x - floor(M*x);
    double yf = M*y - floor(M*y);
    double xc = 0.5 + 0.2*sin(4.615*(x+y));
    double yc = 0.5 + 0.2*sin(4.634*(2*x-3*y));
    double hr = 0.1 + 0.05*sin(22.01*x);
    double hd = hypot(xf-xc, yf-yc);
    if (hd <= hr) { bad = TRUE; }

    /* Add a big value to the slope where we want it to be bad: */
    /* This should set the weight map to zero in those places: */
    if (bad) { dz->c[0] = dz->c[1] = NAN; }
  }

void pst_proc_map_function_20(r2_t p, double *z, r2_t *dz)
  { /* "cpiece" - buried annulus connected by straight ramps: */

    double x = p.c[0], y = p.c[1];

    auto void fix_ramp(double r, double z0, double z1, double t, double w, double d, double s, double *z, r2_t *dz);
      /* Defines the {z,dz} values inside a ramp of width {w} and
        total length {d}. The ramp connects a region at level {z0}
        (bottom) to a region at level {z1} (top), separated by a
        circle with radius {r}. The ramp is located at angle {t} from
        the X-axis. The ramp consists of a flat region with length {s}
        and a slanted plane with length {d-s}. The ramp is directed
        outwards if {d>0}, inwards if {d<0}. The sign of {s} must be
        the same as that of {d}. */

    void fix_ramp(double r, double z0, double z1, double t, double w, double d, double s, double *z, r2_t *dz)
      {
        /* Rotate {x,y} so that the ramp is at +X: */
        double ct = cos(t), st = sin(t);
        double xr = + ct*x + st*y;
        double yr = - st*x + ct*y;

        if ((fabs(yr) <= w/2) && (xr > 0))
          { /* Compute X of beginning of ramp: */
            double xbeg = (d > 0 ? sqrt(r*r - w*w/4) - 0.00001*r : 1.00001*r);
            /* Compute the X position relative to the begin & end of ramp: */
            double u = (xr - xbeg)/d;
            if ((u >= 0) && (u < 1.0))
              { double ubot = s/d;
                if (u < ubot)
                  { if (z != NULL) { (*z) = z0; }
                    if (dz != NULL) { (*dz) = (r2_t){{ 0, 0 }}; }
                  }
                else
                  { double v = (u - ubot)/(1.0 - ubot);
                    if (z != NULL) { (*z) = (1-v)*z0 + v*z1; }
                    if (dz != NULL)
                      { /* Compute gradient relative to {xr,yr}: */
                        double dzdxr = (z1 - z0)/(d - s);
                        /* Rotate gradient to original position: */
                        dz->c[0] = dzdxr*ct;
                        dz->c[1] = dzdxr*st;
                      }
                  }
              }
          }
      }

    double zA = 0.10; /* Background level. */
    double zB = 0.00; /* Annulus level. */
    double zC = 0.05; /* Inner platform level. */

    /* Initialize {*dz} with a plane: */
    if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }

    /* Determine the basic floor level: */
    double rA = 0.75; /* Outer annulus radius. */
    double rB = 0.50; /* Outer annulus radius. */
    double r2 = x*x + y*y;
    if (z != NULL)
      { if (r2 > rA*rA)
          { (*z) = zA; }
        else if (r2 < rB*rB)
          { (*z) = zC; }
        else
          { (*z) = zB; }
      }

    /* Excavate outer ramp: */
    double tA =  0.80; /* Position in radians. */
    double wA =  0.20; /* Width. */
    double dA = +0.30; /* Depth. */
    double sA = +0.05; /* Threshold depth. */
    fix_ramp(rA, zB, zA, tA, wA, dA, sA, z, dz);

    /* Excavate inner ramp: */
    double tB =  2.50; /* Position in radians. */
    double wB =  0.16; /* Width. */
    double dB = -0.24; /* Depth. */
    double sB = -0.04; /* Threshold depth. */
    fix_ramp(rB, zB, zC, tB, wB, dB, sB, z, dz);

  }

void pst_proc_map_function_21(r2_t p, double *z, r2_t *dz)
  { /* "cplat3" - circular platforms with smooth edges. */
    r2_t ctr[3] = { (r2_t){{+0.3,+0.4}}, (r2_t){{-0.5,-0.2}}, (r2_t){{+0.4,-0.5}} }; /* Centers. */
    double R[3] = { 0.5, 0.4, 0.3 };      /* Nominal radii. */
    double HS[3] = { 0.05, 0.02, 0.03 };  /* Half-widthd of shoulders. */
    double zDif[3] = { 0.80, 0.60, 0.70 };   /* Heights. */
    int32_t k;
    if (z != NULL) { (*z) = 0; }
    if (dz != NULL) { (*dz) = (r2_t){{0,0}}; }
    for (k = 0; k < 3; k++)
      {
        r2_t vk; r2_sub(&p, &(ctr[k]), &vk);
        double zk;
        r2_t dzk;
        pst_proc_map_function_round_platform(&vk, R[k], HS[k], &zk, (dz == NULL ? NULL : &dzk));
        if (z != NULL) { (*z) += zDif[k]*zk; }
        if (dz != NULL)
          { dz->c[0] += zDif[k]*dzk.c[0];
            dz->c[1] += zDif[k]*dzk.c[1];
          }
      }
  }

void pst_proc_map_function_22(r2_t p, double *z, r2_t *dz)
  { /* "plat3a" - triangular wall with smooth edges. */
    double R[2] = { 0.9, 0.4 };     /* Nominal radii. */
    double S[2] = { 0.05, 0.03 };   /* Half-widthd of shoulders. */
    double zDif[2] = { +0.80, -0.40 }; /* Heights. */
    double T[2] = { 0, M_PI/3 };    /* Tilts. */
    int32_t k;
    if (z != NULL) { (*z) = 0; }
    if (dz != NULL) { (*dz) = (r2_t){{0,0}}; }
    for (k = 0; k < 2; k++)
      { double zk;
        r2_t dzk;
        pst_proc_map_function_polygonal_platform(&p, 3, R[k], T[k] + M_PI/12, S[k], &zk, (dz == NULL ? NULL : &dzk));
        if (z != NULL) { (*z) += zDif[k]*zk; }
        if (dz != NULL)
          { dz->c[0] += zDif[k]*dzk.c[0];
            dz->c[1] += zDif[k]*dzk.c[1];
          }
      }
  }

void pst_proc_map_function_23(r2_t p, double *z, r2_t *dz)
  { /* "plat5a" - pentagonal platform with smooth edges. */
    double R = 0.6;
    double S = 0.2;
    double T = M_PI/15.0;
    pst_proc_map_function_polygonal_platform(&p, 5, R, T, S, z, dz);
  }

void pst_proc_map_function_24(r2_t p, double *z, r2_t *dz)
  { /* "cplatv" - circular platform with sharp edges. */

    r2_t ctr = (r2_t){{ 0.0, 0.0 }};
    double zMin = 0.05; /* Z-level outside. */
    double zMax = 0.95; /* Z-level inside. */
    double zDif = zMax - zMin;
    double R = 0.8;     /* Nominal radius. */
    double HS = 0.0;    /* Half-widthd of shoulder. */

    r2_t v; r2_sub(&p, &ctr, &v);
    double zt;
    r2_t dzt;
    pst_proc_map_function_round_platform(&v, R, HS, &zt, &dzt);
    if (z != NULL) { (*z) = zMin + zDif*zt; }
    if (dz != NULL)
      { dz->c[0] = zDif*dzt.c[0];
        dz->c[1] = zDif*dzt.c[1];
      }
  }

void pst_proc_map_function_25(r2_t p, double *z, r2_t *dz)
  { /* "cplats" - circular platform with smooth edges. */

    r2_t ctr = (r2_t){{ 0.0, 0.0 }};
    double zMin = 0.05; /* Z-level outside. */
    double zMax = 0.95; /* Z-level inside. */
    double zDif = zMax - zMin;
    double R = 0.70;     /* Nominal radius. */
    double HS = 0.25;   /* Half-width of shoulder. */

    r2_t v; r2_sub(&p, &ctr, &v);
    double zt;
    r2_t dzt;
    pst_proc_map_function_round_platform(&v, R, HS, &zt, &dzt);
    if (z != NULL) { (*z) = zMin + zDif*zt; }
    if (dz != NULL)
      { dz->c[0] = zDif*dzt.c[0];
        dz->c[1] = zDif*dzt.c[1];
      }
  }

void pst_proc_map_function_26(r2_t p, double *z, r2_t *dz)
  { /* "qtbell" - quartic bell. */
    r2_t ctr = (r2_t){{ 0.0, 0.0 }};
    double zMin = 0.05; /* Z-level outside. */
    double zMax = 0.95; /* Z-level inside. */
    double R = 0.8;     /* Nominal radius. */

    double X = (p.c[0] - ctr.c[0])/R; double dXdx = 1/R;
    double Y = (p.c[1] - ctr.c[1])/R; double dYdy = 1/R;
    double zDif = zMax - zMin;

    double r2 = X*X + Y*Y;
    if (r2 >= 1.0)
      { if (z != NULL) { (*z) = zMin; }
        if (dz != NULL) { (*dz) = (r2_t){{0, 0}}; }
      }
    else
      { double dr2dx = 2*X*dXdx; double dr2dy = 2*Y*dYdy;
        double F2 = 1.0 - r2; double dF2dx = - dr2dx; double dF2dy = - dr2dy;
        double F4 = F2*F2; double dF4dx = 2*F2*dF2dx; double dF4dy = 2*F2*dF2dy;
        if (z != NULL) { (*z) = zMin + zDif*F4; }
        if (dz != NULL)
          { dz->c[0] = zDif*dF4dx;
            dz->c[1] = zDif*dF4dy;
          }
      }
  }

void pst_proc_map_function_27(r2_t p, double *z, r2_t *dz)
  { /* "fracto" - fractalish round mountain. */

    r2_t ctr = (r2_t){{ 0.0, 0.0 }};
    double zMin = 0.05; /* Z-level outside. */
    double zMax = 0.95; /* Z-level inside. */
    double R = 0.8;     /* Nominal radius. */

    r2_t v; r2_sub(&p, &ctr, &v); r2_scale(1/R, &v, &v);
    double zt;
    r2_t dzt;
    double eps = 0.0005; /* Enough for 1000x1000 images. Should be a parameter... */
    double seed = (sqrt(5)-1)/2; /* Sort of random... */
    pst_proc_map_fractal_cone(&v, eps, seed, &zt, &dzt);
    double zDif = zMax - zMin;
    if (z != NULL) { (*z) = zMin + zDif*zt; }
    if (dz != NULL)
      { dz->c[0] = zDif*dzt.c[0]/R;
        dz->c[1] = zDif*dzt.c[1]/R;
      }
  }

/* PARAMETRIZED FIELDS */

void pst_proc_map_function_affine(r2_t *p, r2_t *a, double fa, r2_t *b, double fb, r2_t *c, double fc, double *z, r2_t *dz)
  {
    r2_t vca; r2_sub(a, c, &vca);
    r2_t vcb; r2_sub(b, c, &vcb);
    /* Compute the 2x2 matrix {N} that maps {vca} to {(1,0)} to and {vcb} to {(0,1)} : */
    r2x2_t M = (r2x2_t){{{ vca.c[0], vca.c[1] }, {vcb.c[0], vcb.c[1] }}};
    r2x2_t N; r2x2_inv(&M, &N);
    /* Get the barycentric coordinates of {p} relative to {a,b,c}: */
    r2_t vcp; r2_sub(p, c, &vcp);
    r2_t bar; r2x2_map_row(&vcp, &N, &bar);
    double ua = bar.c[0]; double duadx = N.c[0][0], duady = N.c[1][0];
    double ub = bar.c[1]; double dubdx = N.c[0][1], dubdy = N.c[1][1];
    double uc = 1 - ua - ub; double ducdx = - (duadx + dubdx), ducdy = - (duady + dubdy);
    /* Barycentric interpolation of {fa,fb,fc}: */
    if (z != NULL) { (*z) = fa*ua + fb*ub + fc*uc; }
    if (dz != NULL)
      { dz->c[0] = fa*duadx + fb*dubdx + fc*ducdx;
        dz->c[1] = fa*duady + fb*dubdy + fc*ducdy;
      }
  }

void pst_proc_map_function_wave(r2_t *p, r2_t *f, double phase, double *z, r2_t *dz)
  {
    double t = 2*M_PI*r2_dot(p, f) + phase;
    double fm = r2_norm(f);
    double A = 1.5/fmax(0.5, 2.0*M_PI*fm);
    if (z != NULL) { (*z) = A*sin(t); }
    if (dz != NULL)
      { double ct = A*2*M_PI*cos(t);
        dz->c[0] = ct*f->c[0];
        dz->c[1] = ct*f->c[1];
      }
  }

void pst_proc_map_function_babel(r2_t *p, double RI, double RO, double N, double EF, double *z, r2_t *dz)
  { /* Generic tower of Babel: */
    double x = p->c[0], y = p->c[1];

    double zBot = 0.1;  /* Height of ground floor. */
    double zTop = 0.80; /* Height of top platform. */

    double DR = (RO - RI)/(N+1);  /* Width of ramp. */
    double zDif = (zTop - zBot)/N;  /* Height difference between adjacent turns. */
    double EP = EF*DR;            /* Width of shoulder. */

    /* Convert {x,y} to polar {r,t}: */
    double r2 = x*x + y*y;
    double r = sqrt(r2);              /* Distance from center. */
    double t = atan2(y, x)/(2*M_PI);  /* Argument of {(x,y)} in turns. */

    /* Compute height {mz} of ramp, ignoring shoulders and Z-clipping: */
    double rpt = RI + (1 - t)*DR;     /* Radius of platform at this {t}. */
    double s = r - rpt;               /* Radial distance from platform. */
    int32_t k = (int32_t)floor(s/DR);              /* Index of ramp turn (0 = nearest to center). */
    double ph = fmax(0,s - k*DR);     /* Distance from inner edge of ramp, ignoring shoulder. */
    double mz = zTop - zDif*(1 + k - t);  /* Height without shoulders or clipping. */

    if (mz >= zTop)
      { /* Central platform: */
        if (z != NULL) { (*z) = zTop; }
        if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }
      }
    else if ((mz <= zBot) && (ph > EP))
      { /* The ramp here is below ground and we are not on the shoulder, so we are on the ground: */
        if (z != NULL) { (*z) = zBot; }
        if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }
      }
    else if (mz <= zBot - zDif)
      { /* The ramp here is buried so deeply that there is not even the shoulder, so we are on the ground: */
        if (z != NULL) { (*z) = zBot; }
        if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }
      }
    else
      { /* We are on the ramp or shoulder (including the outermost shoulder): */

        double drdx = x/r;
        double drdy = y/r;

        double dtdx = -y/r2/(2*M_PI);
        double dtdy = +x/r2/(2*M_PI);

        double dmdx = zDif*dtdx;
        double dmdy = zDif*dtdy;

        /* First compute {z,dz} ignoring the shoulder fill: */
        if (mz < zBot)
          { /* Ramp is buried so we would be on the ground: */
            if (z != NULL) { (*z) = zBot; }
            if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }
          }
        else
          { /* We would be on the ramp: */
            if (z != NULL) { (*z) = mz; }
            if (dz != NULL)
              { dz->c[0] = dmdx;
                dz->c[1] = dmdy;
              }
          }

        if (ph < EP)
          { /* We are on the shoulder region, add the shoulder fill: */

            /* Get the relative position {w} across the shoulder, in [0_1]: */
            double w = ph/EP;

            /* Compute the the derivatives of {w}: */
            double dsdx = drdx + DR*dtdx;
            double dsdy = drdy + DR*dtdy;

            double dwdx = dsdx/EP;
            double dwdy = dsdy/EP;

            /* We need a cubic Hermite shoulder {h(w)} of unit width and height: */
            double h = (2*w + 1)*(w - 1)*(w - 1);
            double dhdw = 6*w*(w - 1);
            double dhdx = dhdw*dwdx;
            double dhdy = dhdw*dwdy;

            /* Now get the height {a} of the shoulder: */
            double a, dadx, dady;
            if (mz + zDif > zTop)
              { /* Innermost shoulder: */
                a = zTop - mz; dadx = -dmdx; dady = -dmdy;
              }
            else if (mz < zBot)
              { /* Outermost shoulder. */
                a = mz + zDif - zBot; dadx = dmdx; dady = dmdy;
              }
            else
              { /* Full-height shoulder: */
                a = zDif; dadx = dady = 0;
              }

            /* Now add the shoulder of height {a} to {mz}: */
            if (z != NULL) { (*z) += a*h; }
            if (dz != NULL)
              { dz->c[0] += dadx*h + a*dhdx;
                dz->c[1] += dady*h + a*dhdy;
              }
          }
      }
  }

void pst_proc_map_function_cubic_ramp(double x, double *z, double *dz)
  {
    if (z != NULL) { (*z) = 0.25*x*(3 - x*x) + 0.50; }
    if (dz != NULL) { (*dz) = 0.75*(1 - x*x); }
  }

void pst_proc_map_function_round_platform(r2_t *p, double R, double HS, double *z, r2_t *dz)
  {
    double x = p->c[0], y = p->c[1];

    demand(R > 0, "bad radius {R}");
    demand((HS >= 0) && (HS <= R), "bad {HS}");

    /* Default derivative: */
    if (dz != NULL) { (*dz) = (r2_t){{0, 0}}; }

    double RI = R - HS;
    double RO = R + HS;

    double r2 = x*x + y*y;
    if (r2 <= RI*RI)
      { if (z != NULL) { (*z) = 1.0; } }
    else if (r2 >= RO*RO)
      { if (z != NULL) { (*z) = 0.0; } }
    else
      { assert(HS > 0);
        double r = sqrt(r2);
        /* Compute {z} as a function of position within smoothed band: */
        double t = (r - R)/HS; /* Range {[-1 _ +1]}. */
        double h, dh_dt;
        pst_proc_map_function_cubic_ramp(t, &h, (dz == NULL ? NULL : &dh_dt));
        if (z != NULL) { (*z) = 1.0 - h; }
        if (dz != NULL)
          { /* Compute derivative w.r.t. {x,y}: */
            double dh_dr = - dh_dt/HS;
            dz->c[0] = dh_dr*(x/r);
            dz->c[1] = dh_dr*(y/r);
          }
      }
  }

void pst_proc_map_function_polygonal_platform(r2_t *p, int32_t N, double R, double tilt, double S, double *z, r2_t *dz)
  {
    double x = p->c[0], y = p->c[1];

    demand(R > 0, "bad radius {R}");
    demand(S >= 0, "bad {S}");
    demand(N >= 3, "bad {N}");

    /* Default derivative: */
    if (dz != NULL) { (*dz) = (r2_t){{0, 0}}; }

    double RI = R*cos(M_PI/N); /* Inscribed circle radius. */
    double RO = R + S;  /* Outer radius including smooth transition. */

    double r2 = x*x + y*y;
    if (r2 >= RO*RO)
      { if (z != NULL) { (*z) = 0; } }
    else if (r2 <= RI*RI)
      { if (z != NULL) { (*z) = 1; } }
    else
      {
        /* Compute sectors to turn: */
        double ang = atan2(y,x) - tilt; /* Angle relative to ref vertex. */
        int32_t nw = (int32_t)floor(N*ang/(2*M_PI));  /* Sectors to turn. */

        /* Convert to a coordinate system {u,v} where the mid-side is at {v=0}: */
        double hw = M_PI/N;               /* Half angular width of each sector. */
        double t = (1 + 2*nw)*hw + tilt;  /* Total angle to turn clockwise (in radians). */
        double ct = cos(-t);
        double st = sin(-t);

        double u = + x*ct - y*st;
        double v = + x*st + y*ct;

        /* Fold negative {v}: */
        double sv = (v < 0 ? -1 : +1);
        v = sv*v;

        double HL = R*sin(hw); /* Half-side of outline. */

        if (u <= 0)
          { fprintf(stderr, "p = ( %+8.5f %+8.5f )  ang = %+8.2f  nw = %+2d", x, y, ang*180/M_PI, nw);
            fprintf(stderr, "  t = %+8.2f  uv = ( %+8.5f %+8.5f )  auv = %+8.2f\n", t*180/M_PI, u, v, atan2(v,u)*180/M_PI);
          }

        assert(u > 0);
        assert(v/u <= 1.00001*HL/RI);

        if (u <= RI)
          { if (z != NULL) { (*z) = 1; } }
        else if (u >= RI + S)
          { if (z != NULL) { (*z) = 0; } }
        else
          { if (v <= HL)
              { /* Ramp on straight part of outline: */
                double m = 2*(u - RI)/S - 1; /* Position along ramp, in {[-1_+1]}. */
                double h, dh_dm;
                pst_proc_map_function_cubic_ramp(m, &h, (dz == NULL ? NULL : &dh_dm));
                if (z != NULL) { (*z) = 1.0 - h; }
                if (dz != NULL)
                  { /* Compute derivative w.r.t. {x,y}: */
                    double dz_du = - 2*dh_dm/S;
                    dz->c[0] = + dz_du*ct;
                    dz->c[1] = - dz_du*st;
                  }
              }
            else
              { /* Ramp on corner: */
                double q = hypot(u - RI, v - HL);
                if (q >= S)
                  { if (z != NULL) { (*z) = 0; } }
                else
                  { double m = 2*q/S - 1; /* Position along ramp, in {[-1_+1]}. */
                    double h, dh_dm;
                    pst_proc_map_function_cubic_ramp(m, &h, (dz == NULL ? NULL : &dh_dm));
                    (*z) = 1.0 - h;
                    if (dz != NULL)
                      { /* Compute derivative w.r.t. {x,y}: */
                        double dz_dq = - 2*dh_dm/S;
                        double dz_du = dz_dq*(u - RI)/q;
                        double dz_dv = dz_dq*(v - HL)/q;
                        dz->c[0] = + dz_du*ct + sv*dz_dv*st;
                        dz->c[1] = - dz_du*st + sv*dz_dv*ct;
                      }
                  }
              }
          }
      }
  }
