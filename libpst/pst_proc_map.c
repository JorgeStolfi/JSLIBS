/* See pst_proc_map.h */
/* Last edited on 2025-01-20 08:02:16 by stolfi */

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

/* IMPLEMENTATIONS */
  
void pst_proc_map_compute_height
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    uint32_t NS,
    double ws[],
    double xyScale,
    double maxGrad,
    double *z_P,
    double *w_P
  )
  { 
    demand(isfinite(xyScale) && (xyScale > 0), "invalid {xyScale}");
    demand(isfinite(maxGrad) && (maxGrad > 0), "invalid maxGrad");
    demand(((NS%2) == 1) && (NS >= 3), "sampoint count {NS} must be odd and at least 3");

    int32_t HS = (int32_t)NS/2;
    double step = 1.0/(HS+1.0)/xyScale; /* Displacement between sampoiints. */
    
    /* Accumulators for average height: */
    double sum_wkzk = 0;
    double sum_wk = 0;
    double zprev[NS]; /* Height value in previous sampoint row. */
    
    double max_rdz = 0; /* Max relative gradient seen. */
    
    for (int32_t ky = -HS; ky <= +HS; ky++)
      { for (int32_t kx = -HS; kx <= +HS; kx++)
          { /* Generate a sampoint {(xk,yk)} around {p}: */
            double xk = p.c[0] + ((double)kx)*step;
            double yk = p.c[1] + ((double)ky)*step;
            r2_t pk = (r2_t){{ xk, yk }};
            /* Get the sampoint weight {wk}: */
            double wk = ws[HS+kx]*ws[HS+ky];
            /* Get height value, ignore gradient: */
            double zk;
            func(pk, &zk, NULL);
            if (! isfinite(zk)) { max_rdz = +INF; zk = 0.0; }
            /* Accumulate height value: */
            sum_wkzk += wk*zk;
            sum_wk += wk;
            /* Check gradient, for weight's sake: */
            double dzx = 0, dzy = 0;
            if (kx > 0) 
              { /* Check {X} derivative: */
                dzx = fabs(zk - zprev[HS+kx-1])/step;
              }
            if (ky > 0)
              { /* Check {Y} derivative: */
                dzy = fabs(zk - zprev[HS+kx])/step;
              }
            zprev[HS+kx] = zk;
            double rdz = hypot(dzx,dzy)/maxGrad;
            max_rdz = fmax(max_rdz, rdz);
          }
      }
    /* Return average {z} and weight: */
    (*z_P) = (float)(sum_wk == 0 ? 0.0: sum_wkzk/sum_wk);
    (*w_P) = (float)(max_rdz >= 1.0 ? 0.0 : 1 - max_rdz*max_rdz);
  }  

void pst_proc_map_compute_gradient
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    uint32_t NS,
    double ws[],
    double xyScale,
    bool_t numGrad,
    double maxGrad,
    double maxGDiff,
    r2_t *dz_P,
    float *w_P
  )
  {
    demand(isfinite(xyScale) && (xyScale > 0), "invalid {xyScale}");
    demand(isfinite(maxGrad) && (maxGrad > 0), "invalid maxGrad");
    demand(isfinite(maxGDiff) && (maxGDiff > 0), "invalid maxGDiff");
    demand(((NS%2) == 1) && (NS >= 3), "sampoint count {NS} must be odd and at least 3");

    int32_t HS = (int32_t)NS/2;
    double step = 1.0/(HS+1.0)/xyScale; /* Displacement between sampoiints. */
    
    /* Accumulators for average analytic gradient: */
    r2_t sum_wa_dz = (r2_t){{ 0, 0 }};
    double sum_wa = 0;
    
    /* Accumulators for average numeric gradient: */
    r2_t sum_wn_dz = (r2_t){{ 0, 0 }};
    r2_t sum_wn = (r2_t){{ 0, 0 }};
    double zprev[NS]; /* Height values on previous sampoint row. */
     
    double max_rdz = 0; /* Max relative gradient seen. */
    
    for (int32_t ky = -HS; ky <= +HS; ky++)
      { for (int32_t kx = -HS; kx <= +HS; kx++)
          { /* Generate a sampoint {(xk,yk)} around {p}: */
            double xk = p.c[0] + ((double)kx)*step;
            double yk = p.c[1] + ((double)ky)*step;
            r2_t pk = (r2_t){{ xk, yk }};
            /* Get the sampoint weight {wk}: */
            double wk = ws[HS+kx]*ws[HS+ky];
            /* Get height value and analytic derivatives: */
            double zk; r2_t dazk;
            func(pk, &zk, &dazk);
            if (! isfinite(zk)) { max_rdz = +INF; zk = 0.0; }
            if (! isfinite(dazk.c[0])) { max_rdz = +INF; dazk.c[0] = 0.0; }
            if (! isfinite(dazk.c[1])) { max_rdz = +INF; dazk.c[1] = 0.0; }
            /* Accumulate the analytic gradient: */
            sum_wa_dz.c[0] += wk*dazk.c[0];
            sum_wa_dz.c[1] += wk*dazk.c[1];
            sum_wa += wk;
            double dnzxk = 0, dnzyk = 0;
            if (kx > 0)
              { /* Compute and accumulate numeric {X} derivative: */
                dnzxk = (zk - zprev[HS+kx-1])/step;
                double wnx = (ws[HS+kx] + ws[HS+kx-1])/2;
                sum_wn_dz.c[0] += wnx*dnzxk;
                sum_wn.c[0] += wnx;
              }
            if (ky > 0)
              { /* Compute and accumulate numeric {Y} derivative: */
                dnzyk = (zk - zprev[HS+kx])/step;
                double wny = (ws[HS+ky] + ws[HS+ky-1])/2;
                sum_wn_dz.c[1] += wny*dnzyk;
                sum_wn.c[1] += wny;
              }
            zprev[HS+kx] = zk;
            /* Update max relative gradient: */
            /* If {zk} was invalid this is nonsense but does not matter: */
            double rdz = hypot(dnzxk, dnzyk)/maxGrad;
            max_rdz = fmax(max_rdz, rdz);
         }
     }
        
    /* Compute mean numeric gradient: */
    double dnzx_avg = (sum_wn.c[0] == 0 ? 0.0 : sum_wn_dz.c[0]/sum_wn.c[0]);
    double dnzy_avg = (sum_wn.c[1] == 0 ? 0.0 : sum_wn_dz.c[1]/sum_wn.c[1]);
    r2_t dnz_avg = (r2_t){{ dnzx_avg, dnzy_avg }};

    double dazx_avg = (sum_wa == 0 ? 0.0: sum_wa_dz.c[0]/sum_wa);
    double dazy_avg = (sum_wa == 0 ? 0.0: sum_wa_dz.c[1]/sum_wa);
    r2_t daz_avg = (r2_t){{ dazx_avg, dazy_avg }};
        
    (*dz_P) = (numGrad ? dnz_avg : daz_avg);
            
    /* Compute weight from {max_rdz} and num/analytic discrepancy: */
    double w_avg;
    if (max_rdz >= 1.0)
      { w_avg = 0.0; }
    else
      { double rg = hypot(dnzx_avg  - dazx_avg, dnzy_avg - dazy_avg)/maxGDiff;
        w_avg = (rg >= 1.0 ? 0.0 : (1 - rg*rg));
      }
    (*w_P) = (float)w_avg;
  }

float_image_t* pst_proc_map_make_height_map
  ( pst_proc_map_zfunc_t *func,
    int32_t NX,
    int32_t NY,
    uint32_t NS,
    double ws[],
    double maxGrad
  )
  {
    int32_t NC = 2;
    float_image_t *IZ = float_image_new(NC, NX+1, NY+1);

    int32_t NMIN = (NX < NY ? NX : NY);  /* Smallest dimension of image. */
    double xyScale = ((double)NMIN)/2.0;  /* Number of grid pixels per {func} domain unit. */
    /* Center of image domain (in pixels, from bottom left corner): */
    double cX = 0.5*NX;
    double cY = 0.5*NY;
    
    /* Compute height at each grid corner: */
    for (int32_t y = 0; y <= NY; y++)
      { for (int32_t x = 0; x <= NX; x++)
          { /* Compute function domain coordinates {(xp,yp)} of grid corner: */
            double xp = (x - cX)/xyScale;
            double yp = (y - cY)/xyScale;
            r2_t p = (r2_t){{ xp, yp }};

            /* Compute height at grid corner: */
            double z, w;
            pst_proc_map_compute_height(func, p, NS, ws, xyScale, maxGrad, &z, &w);

            /* Scale height from {func} units to pixels: */
            float_image_set_sample(IZ, 0, x, y, (float)(z*xyScale));
            float_image_set_sample(IZ, 1, x, y, (float)w);
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
          { /* Compute function domain coordinates {(xp,yp)} of pixel center: */
            double xp = (x + 0.5 - cX)/xyScale;
            double yp = (y + 0.5 - cY)/xyScale;
            r2_t p = (r2_t){{ xp, yp }};

            /* Compute gradient at pixel center: */
            r2_t dz; float w;
            pst_proc_map_compute_gradient
              ( func, p, NS, ws, xyScale, numGrad, maxGrad, maxGDiff, &dz, &w );
            float_image_set_sample(IG, 0, x, y, (float)dz.c[0]); 
            float_image_set_sample(IG, 1, x, y, (float)dz.c[1]);
            float_image_set_sample(IG, 2, x, y, (float)w);
          }
      }
    return IG;
  }

pst_proc_map_zfunc_t *pst_proc_map_function_generic(int32_t n)
  {
    switch(n)
      { case  0: return &pst_proc_map_function_00; /* "zeflat" Constant function (zero gradient). */
        case  1: return &pst_proc_map_function_01; /* "ramp10" Linear ramp in the X direction. */
        case  2: return &pst_proc_map_function_02; /* "ramp01" Linear ramp in the Y direction. */
        case  3: return &pst_proc_map_function_03; /* "ramp11" Linear ramp in the XY diagonal direction. */
        case  4: return &pst_proc_map_function_04; /* "parabo" Slanted elliptical parabolic hump. */
        case  5: return &pst_proc_map_function_05; /* "spdom1" Buried sphere. */
        case  6: return &pst_proc_map_function_06; /* "pyram5" Pentagonal pyramid. */
        case  7: return &pst_proc_map_function_07; /* "rtcone" Conical mound. */
        case  8: return &pst_proc_map_function_08; /* "ripple" Circular waves. */
        case  9: return &pst_proc_map_function_09; /* "sbabel" Tower of Babel with smooth shoulders. */
        case 10: return &pst_proc_map_function_10; /* "hcliff" Diverging parabolic ramps with cliff. */
        case 11: return &pst_proc_map_function_11; /* "sinw01" Sinusoidal wave, low freq. */
        case 12: return &pst_proc_map_function_12; /* "cosw01" Co-sinusoidal wave, low freq. */
        case 13: return &pst_proc_map_function_13; /* "sinw02" Sinusoidal wave, med freq. */
        case 14: return &pst_proc_map_function_14; /* "cosw02" Co-sinusoidal wave, med freq. */
        case 15: return &pst_proc_map_function_15; /* "sinw03" Sinusoidal wave, high freq. */
        case 16: return &pst_proc_map_function_16; /* "cosw03" Co-sinusoidal wave, high freq. */
        case 17: return &pst_proc_map_function_17; /* "cbabel" Tower of Babel with cliff. */
        case 18: return &pst_proc_map_function_18; /* "cbramp" Cubic ramp with cliff on three sides. */
        case 19: return &pst_proc_map_function_19; /* "holes1" Buried sphere with many holes in slope map. */
        case 20: return &pst_proc_map_function_20; /* "cpiece" Three nested stages connected by ramps. */
        case 21: return &pst_proc_map_function_21; /* "cplat3" Three flat round platforms with smooth edges. */
        case 22: return &pst_proc_map_function_22; /* "plat3a" Triangular wall with smooth edges. */
        case 23: return &pst_proc_map_function_23; /* "plat5a" Pentagonal platform with smooth edges. */
        case 24: return &pst_proc_map_function_24; /* "cplatv" Circular platform with sharp edges. */
        case 25: return &pst_proc_map_function_25; /* "cplats" Circular platform with smooth edges. */
        case 26: return &pst_proc_map_function_26; /* "qtbell" Quartic bell. */
        case 27: return &pst_proc_map_function_27; /* "fracto" Fractalish round montain. */
        default:
          fprintf(stderr, "bad {n} %d\n", n);
          assert(FALSE);
      }
  }

char* pst_proc_map_function_generic_name(int32_t n)
  {
    switch(n)
      { case  0: return "zeflat"; /* Constant function (zero gradient). */
        case  1: return "ramp10"; /* Linear ramp in the X direction. */
        case  2: return "ramp01"; /* Linear ramp in the Y direction. */
        case  3: return "ramp11"; /* Linear ramp in the XY diagonal direction. */
        case  4: return "parabo"; /* Slanted elliptical parabolic hump. */
        case  5: return "spdom1"; /* Buried sphere. */
        case  6: return "pyram5"; /* Pentagonal pyramid. */
        case  7: return "rtcone"; /* Conical mound. */
        case  8: return "ripple"; /* Circular waves. */
        case  9: return "sbabel"; /* Tower of Babel with smooth shoulders. */
        case 10: return "hcliff"; /* Diverging parabolic ramps with cliff. */
        case 11: return "sinw01"; /* Sinusoidal wave, low freq. */
        case 12: return "cosw01"; /* Co-sinusoidal wave, low freq. */
        case 13: return "sinw02"; /* Sinusoidal wave, med freq. */
        case 14: return "cosw02"; /* Co-sinusoidal wave, med freq. */
        case 15: return "sinw03"; /* Sinusoidal wave, high freq. */
        case 16: return "cosw03"; /* Co-sinusoidal wave, high freq. */
        case 17: return "cbabel"; /* Tower of Babel with cliff. */
        case 18: return "cbramp"; /* Cubic ramp with cliff on three sides. */
        case 19: return "holes1"; /* Buried sphere with many holes in slope map. */
        case 20: return "cpiece"; /* Three nested stages connected by ramps. */
        case 21: return "cplat3"; /* Three flat round platforms with smooth edges. */
        case 22: return "plat3a"; /* Triangular wall with smooth edges. */
        case 23: return "plat5a"; /* Pentagonal platform with smooth edges. */
        case 24: return "cplatv"; /* Circular platform with sharp edges. */
        case 25: return "cplats"; /* Circular platform with smooth edges. */
        case 26: return "qtbell"; /* Quartic bell. */
        case 27: return "fracto"; /* Fractalish round montain. */
        default:
          fprintf(stderr, "bad {n} %d\n", n);
          assert(FALSE);
      }
  }

void pst_proc_map_function_00(r2_t p, double *z, r2_t *dz)
  { /* "zeflat" - constant function: */
    double zval = 0.50;
    (*z) = zval;
    if (dz != NULL) { dz->c[0] = dz->c[1] = 0; }
  }

void pst_proc_map_function_01(r2_t p, double *z, r2_t *dz)
  { /* "ramp10" - linear ramp along X: */
    double x = p.c[0], y = p.c[1];
    double Cx = 0.5, Cy = 0.0;
    (*z) = Cx*x + Cy*y;
    if (dz != NULL) { dz->c[0] = Cx; dz->c[1] = Cy; }
  }

void pst_proc_map_function_02(r2_t p, double *z, r2_t *dz)
  { /* "ramp01" - linear ramp along Y: */
    double x = p.c[0], y = p.c[1];
    double Cx = 0.0, Cy = 0.5;
    (*z) = Cx*x + Cy*y;
    if (dz != NULL) { dz->c[0] = Cx; dz->c[1] = Cy; }
  }

void pst_proc_map_function_03(r2_t p, double *z, r2_t *dz)
  { /* "ramp11" - linear ramp along diagonal: */
    double x = p.c[0], y = p.c[1];
    double Cx = 1/4.0, Cy = 1/3.0;
    (*z) = Cx*x + Cy*y;
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
    (*z) = 1.0 - A*S*S - B*T*T;
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
    (*z) = zMin;
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
            (*z) = F;
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
    (*z) = zMin;
    if (dz != NULL){ dz->c[0] = dz->c[1] = 0.00; }
    /* Where the pyramid is above the ground, set {*z,*dz} to it: */
    double R = 0.75; /* Pyramid radius at z = 0. */
    double H = 0.80; /* Pyramid height. */
    double D2 = x*x + y*y;
    if (D2 == 0)
      { /* At the tip, assume slope = 0: */
        (*z) = H;
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
                (*z) = F;
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
    (*z) = zMin;
    if (dz != NULL){ dz->c[0] = dz->c[1] = 0.00; }
    double R = 0.75; /* Cone radius at z = 0. */
    double H = 0.80; /* Cone height. */
    double D2 = x*x + y*y;
    if (D2 == 0)
      { /* At the tip, assume slope = 0: */
        (*z) = H;
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
                (*z) = F;
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
    (*z) = F;
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
        (*z) = 0;
        if (dz != NULL)
          { dz->c[0] = 0;
            dz->c[1] = 0;
          }
      }
    else if (T > 0)
      { /* Rising parabolic ramp: */
        (*z) = + A*S*S;
        if (dz != NULL)
          { dz->c[0] = + 2*A*S*Cx;
            dz->c[1] = + 2*A*S*Cy;
          }
      }
    else
      { /* Descending parabolic ramp: */
        (*z) = - A*S*S;
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
        (*z) = 0;
        if (dz != NULL)
          { dz->c[0] = 0;
            dz->c[1] = 0;
          }
      }
    else if (u > +1.00)
      { /* Upper floor: */
        (*z) = H;
        if (dz != NULL)
          { dz->c[0] = 0;
            dz->c[1] = 0;
          }
      }
    else
      { /* Ramp: */
        (*z) = 0.25*H*(2 + 3*u - u*u*u);
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

    /* Now add a large value (20.0) to the slopes at random places: */
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
    double Big = 20; /* Hopefully bigger than any actual slope. */
    if (bad) { dz->c[0] += Big; dz->c[1] += Big; }
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
                  { (*z) = z0;
                    if (dz != NULL) { (*dz) = (r2_t){{ 0, 0 }}; }
                  }
                else
                  { double v = (u - ubot)/(1.0 - ubot);
                    (*z) = (1-v)*z0 + v*z1;
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
    if (r2 > rA*rA)
      { (*z) = zA; }
    else if (r2 < rB*rB)
      { (*z) = zC; }
    else
      { (*z) = zB; }

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
    (*z) = 0;
    if (dz != NULL) { (*dz) = (r2_t){{0,0}}; }
    for (k = 0; k < 3; k++)
      {
        r2_t vk; r2_sub(&p, &(ctr[k]), &vk);
        double zk;
        r2_t dzk;
        pst_proc_map_function_round_platform(&vk, R[k], HS[k], &zk, (dz == NULL ? NULL : &dzk));
        (*z) += zDif[k]*zk;
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
    (*z) = 0;
    if (dz != NULL) { (*dz) = (r2_t){{0,0}}; }
    for (k = 0; k < 2; k++)
      { double zk;
        r2_t dzk;
        pst_proc_map_function_polygonal_platform(&p, 3, R[k], T[k] + M_PI/12, S[k], &zk, (dz == NULL ? NULL : &dzk));
        (*z) += zDif[k]*zk;
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
    (*z) = zMin + zDif*zt;
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
    (*z) = zMin + zDif*zt;
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
      { (*z) = zMin;
        if (dz != NULL) { (*dz) = (r2_t){{0, 0}}; }
      }
    else
      { double dr2dx = 2*X*dXdx; double dr2dy = 2*Y*dYdy;
        double F2 = 1.0 - r2; double dF2dx = - dr2dx; double dF2dy = - dr2dy;
        double F4 = F2*F2; double dF4dx = 2*F2*dF2dx; double dF4dy = 2*F2*dF2dy;
        (*z) = zMin + zDif*F4;
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
    (*z) = zMin + zDif*zt;
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
    (*z) = fa*ua + fb*ub + fc*uc;
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
    (*z) = A*sin(t);
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
        (*z) = zTop;
        if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }
      }
    else if ((mz <= zBot) && (ph > EP))
      { /* The ramp here is below ground and we are not on the shoulder, so we are on the ground: */
        (*z) = zBot;
        if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }
      }
    else if (mz <= zBot - zDif)
      { /* The ramp here is buried so deeply that there is not even the shoulder, so we are on the ground: */
        (*z) = zBot;
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
            (*z) = zBot;
            if (dz != NULL) { dz->c[0] = dz->c[1] = 0.00; }
          }
        else
          { /* We would be on the ramp: */
            (*z) = mz;
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
            (*z) += a*h;
            if (dz != NULL)
              { dz->c[0] += dadx*h + a*dhdx;
                dz->c[1] += dady*h + a*dhdy;
              }
          }
      }
  }

void pst_proc_map_function_cubic_ramp(double x, double *z, double *dz)
  {
    (*z) = 0.25*x*(3 - x*x) + 0.50;
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
      { (*z) = 1.0; }
    else if (r2 >= RO*RO)
      { (*z) = 0.0; }
    else
      { assert(HS > 0);
        double r = sqrt(r2);
        /* Compute {z} as a function of position within smoothed band: */
        double t = (r - R)/HS; /* Range {[-1 _ +1]}. */
        double h, dh_dt;
        pst_proc_map_function_cubic_ramp(t, &h, (dz == NULL ? NULL : &dh_dt));
        (*z) = 1.0 - h;
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
      { (*z) = 0; }
    else if (r2 <= RI*RI)
      { (*z) = 1; }
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
          { (*z) = 1; }
        else if (u >= RI + S)
          { (*z) = 0; }
        else
          { if (v <= HL)
              { /* Ramp on straight part of outline: */
                double m = 2*(u - RI)/S - 1; /* Position along ramp, in {[-1_+1]}. */
                double h, dh_dm;
                pst_proc_map_function_cubic_ramp(m, &h, (dz == NULL ? NULL : &dh_dm));
                (*z) = 1.0 - h;
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
                  { (*z) = 0; }
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
