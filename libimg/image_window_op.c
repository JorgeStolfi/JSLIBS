/* See {image_window_op.h}. */
/* Last edited on 2023-09-24 00:30:16 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <image_window_op.h>

void image_window_op_get_window_size
  ( image_window_op_t op, 
    bool_t smoothed,
    int32_t *nxP, 
    int32_t *nyP, 
    int32_t *ictrP
  )
  {
    int32_t nx = 0; /* Window width */
    int32_t ny = 0; /* Window height */
    if (smoothed)
      { nx = 3; ny = 3; }
    else
      { 
        switch(op)
          { 
            case image_window_op_AVERAGE:    nx = 3; ny = 3; break;
            case image_window_op_VARIANCE:   nx = 3; ny = 3; break;
            case image_window_op_DEVIATION:  nx = 3; ny = 3; break;
            case image_window_op_F:          nx = 1; ny = 1; break;
            case image_window_op_FX:         nx = 3; ny = 1; break;
            case image_window_op_FY:         nx = 1; ny = 3; break;
            case image_window_op_FXX:        nx = 3; ny = 1; break;
            case image_window_op_FXY:        nx = 3; ny = 3; break;
            case image_window_op_FYY:        nx = 1; ny = 3; break;
            case image_window_op_FXXY:       nx = 3; ny = 3; break;
            case image_window_op_FXYY:       nx = 3; ny = 3; break;
            case image_window_op_FXXYY:      nx = 3; ny = 3; break;
            case image_window_op_GRADIENT:   nx = 3; ny = 3; break;
            case image_window_op_LAPLACIAN:  nx = 3; ny = 3; break;
            case image_window_op_ORTHICITY:  nx = 3; ny = 3; break;
            case image_window_op_ELONGATION: nx = 3; ny = 3; break;
            case image_window_op_LINF:       nx = 3; ny = 3; break;
            case image_window_op_LINFX:      nx = 3; ny = 3; break;
            case image_window_op_LINFY:      nx = 3; ny = 3; break;
            case image_window_op_LINVAR:     nx = 3; ny = 3; break;
            case image_window_op_LINDEV:     nx = 3; ny = 3; break;
            default: demand(FALSE, "invalid {op}"); 
          }
      }
    (*nxP) = nx;
    (*nyP) = ny;
    assert(nx % 2 == 1);
    assert(ny % 2 == 1);
    (*ictrP) = (ny/2)*nx + (nx/2);
  }

void image_window_op_get_range(image_window_op_t op, bool_t smoothed, bool_t squared, double *loP, double *hiP)
  { 
    /* Determine the natural range {[lo_hi]} of the operator, without squaring: */
    double lo = NAN, hi = NAN;
    if (! squared)
      {
        switch(op)
          {
            case image_window_op_AVERAGE:    lo = 00.0; hi = +1.0;       break;
            case image_window_op_VARIANCE:   lo = 00.0; hi = +0.25;      break;
            case image_window_op_DEVIATION:  lo = 00.0; hi = +0.5;       break;
            case image_window_op_F:          lo = 00.0; hi = +1.0;       break;
            case image_window_op_FX:         lo = -0.5; hi = +0.5;       break;
            case image_window_op_FY:         lo = -0.5; hi = +0.5;       break;
            case image_window_op_FXX:        lo = -2.0; hi = +2.0;       break;
            case image_window_op_FXY:        lo = -0.5; hi = +0.5;       break;
            case image_window_op_FYY:        lo = -2.0; hi = +2.0;       break;
            case image_window_op_FXXY:       lo = -2.0; hi = +2.0;       break;
            case image_window_op_FXYY:       lo = -2.0; hi = +2.0;       break;
            case image_window_op_FXXYY:      lo = -8.0; hi = +8.0;       break;
            case image_window_op_GRADIENT:   
              if (! smoothed)
                { lo = 00.0; hi = +sqrt(0.5); }
              else
                { lo = 00.00; hi = +sqrt(5)/4; }
              break;
            case image_window_op_LAPLACIAN:
              if (! smoothed)
                { lo = -4.0; hi = +4.0; }
              else
                { lo = -2.0; hi = +2.0; } 
              break;
            case image_window_op_ORTHICITY:  lo = -2.0; hi = +2.0;       break;
            case image_window_op_ELONGATION: lo = 00.0; hi = +sqrt(5.0); break;
            case image_window_op_LINF:       lo = 00.0; hi = +1.0;       break;
            case image_window_op_LINFX:      lo = -0.5; hi = +0.5;       break;
            case image_window_op_LINFY:      lo = -0.5; hi = +0.5;       break;
            case image_window_op_LINVAR:     lo = 00.0; hi = +0.25;      break;
            case image_window_op_LINDEV:     lo = 00.0; hi = +0.5;       break;
            default: demand(FALSE, "invalid {op}"); 
          }
      }
    else
      { lo = 0.0;
        switch(op)
          {
            case image_window_op_AVERAGE:    hi = +1.0;     break;
            case image_window_op_VARIANCE:   hi = +0.0625;  break;
            case image_window_op_DEVIATION:  hi = +0.25;    break;
            case image_window_op_F:          hi = +1.0;     break;
            case image_window_op_FX:         hi = +0.25;    break;
            case image_window_op_FY:         hi = +0.25;    break;
            case image_window_op_FXX:        hi = +4.0;     break;
            case image_window_op_FXY:        hi = +0.25;    break;
            case image_window_op_FYY:        hi = +4.0;     break;
            case image_window_op_FXXY:       hi = +4.0;     break;
            case image_window_op_FXYY:       hi = +4.0;     break;
            case image_window_op_FXXYY:      hi = +64.0;    break;
            case image_window_op_GRADIENT:   
              if (! smoothed)
                { hi = +0.5; }
              else
                { hi = +5.0/16; }
              break;
            case image_window_op_LAPLACIAN: 
              if (! smoothed)
                { hi = +16.0; }
              else
                { hi = +4.0; } 
              break;
            case image_window_op_ORTHICITY:  hi = +4.0;     break;
            case image_window_op_ELONGATION: hi = +5.0;     break;
            case image_window_op_LINF:       hi = +1.0;     break;
            case image_window_op_LINFX:      hi = +0.25;    break;
            case image_window_op_LINFY:      hi = +0.25;    break;
            case image_window_op_LINVAR:     hi = +0.0625;  break;
            case image_window_op_LINDEV:     hi = +0.25;    break;
            default: demand(FALSE, "invalid {op}"); 
          }
      }
    assert(! isnan(lo));
    assert(! isnan(hi));
    /* Return results: */
    (*loP) = lo;
    (*hiP) = hi;
  }

double image_window_op_apply
  ( image_window_op_t op, 
    bool_t smoothed, 
    bool_t squared, 
    int32_t ictr, 
    int32_t nx, 
    double smp[]
  )
  {
    double res = 0;
    switch(op)
      { case image_window_op_AVERAGE:
          res = image_window_op_average(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_VARIANCE:
          res = image_window_op_variance(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_DEVIATION:
          if (squared) 
            { res = image_window_op_variance(ictr, nx, smp); }
          else
            { res = image_window_op_deviation(ictr, nx, smp); }
          break;
        case image_window_op_F:
          res = smp[ictr]; /* For both smoothed and unsmoothed. */
          if (squared) { res = res*res; }
          break;
        case image_window_op_FX:
          if (smoothed)
            { res = image_window_op_fx_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_fx_unsmoothed(ictr, nx, smp); }
          if (squared) { res = res*res; }
          break;
        case image_window_op_FY:
          if (smoothed)
            { res = image_window_op_fy_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_fy_unsmoothed(ictr, nx, smp); }
          if (squared) { res = res*res; }
          break;
        case image_window_op_FXX:
          if (smoothed)
            { res = image_window_op_fxx_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_fxx_unsmoothed(ictr, nx, smp); }
          if (squared) { res = res*res; }
          break;
        case image_window_op_FXY:
          res = image_window_op_fxy(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_FYY:
          if (smoothed)
            { res = image_window_op_fyy_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_fyy_unsmoothed(ictr, nx, smp); }
          if (squared) { res = res*res; }
          break;
        case image_window_op_FXXY:
          res = image_window_op_fxxy(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_FXYY:
          res = image_window_op_fxyy(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_FXXYY:
          res = image_window_op_fxxyy(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_GRADIENT:
          if (smoothed)
            { if (squared) 
                { res = image_window_op_gradient_smoothed_squared(ictr, nx, smp); }
              else
                { res = image_window_op_gradient_smoothed(ictr, nx, smp); }
            }
          else
            { if (squared) 
                { res = image_window_op_gradient_unsmoothed_squared(ictr, nx, smp); }
              else
                { res = image_window_op_gradient_unsmoothed(ictr, nx, smp); }
            }
          break;
        case image_window_op_LAPLACIAN:
          if (smoothed)
            { res = image_window_op_laplacian_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_laplacian_unsmoothed(ictr, nx, smp); }
          if (squared) { res = res*res; }
          break;
        case image_window_op_ORTHICITY:
          res = image_window_op_orthicity(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_ELONGATION:
          if (squared) 
            { res = image_window_op_elongation_squared(ictr, nx, smp); }
          else
            { res = image_window_op_elongation(ictr, nx, smp); }
          break;
        case image_window_op_LINF:
          res = image_window_op_average(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_LINFX:
          res = image_window_op_fx_smoothed(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_LINFY:
          res = image_window_op_fy_smoothed(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_LINVAR:
          res = image_window_op_linvar(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_LINDEV:
          if (squared) 
            { res = image_window_op_linvar(ictr, nx, smp); }
          else
            { res = image_window_op_lindev(ictr, nx, smp); }
          break;
        default: demand(FALSE, "invalid {op}"); 
      }
    return res;
  }

image_window_op_t image_window_op_from_string(const char *chop)
  {
    if (strcmp(chop, "average") == 0)
      { return image_window_op_AVERAGE; }
    else if (strcmp(chop, "variance") == 0)
      { return image_window_op_VARIANCE; }
    else if (strcmp(chop, "deviation") == 0)
      { return image_window_op_DEVIATION; }
    else if (strcmp(chop, "f") == 0)
      { return image_window_op_F; }
    else if (strcmp(chop, "fx") == 0)
      { return image_window_op_FX; }
    else if (strcmp(chop, "fy") == 0)
      { return image_window_op_FY; }
    else if (strcmp(chop, "fxx") == 0)
      { return image_window_op_FXX; }
    else if (strcmp(chop, "fxy") == 0)
      { return image_window_op_FXY; }
    else if (strcmp(chop, "fyy") == 0)
      { return image_window_op_FYY; }
    else if (strcmp(chop, "fxxy") == 0)
      { return image_window_op_FXXY; }
    else if (strcmp(chop, "fxyy") == 0)
      { return image_window_op_FXYY; }
    else if (strcmp(chop, "fxxyy") == 0)
      { return image_window_op_FXXYY; }
    else if (strcmp(chop, "gradient") == 0)
      { return image_window_op_GRADIENT; }
    else if (strcmp(chop, "laplacian") == 0)
      { return image_window_op_LAPLACIAN; }
    else if (strcmp(chop, "orthicity") == 0)
      { return image_window_op_ORTHICITY; }
    else if (strcmp(chop, "elongation") == 0)
      { return image_window_op_ELONGATION; }
    else if (strcmp(chop, "linf") == 0)
      { return image_window_op_LINF; }
    else if (strcmp(chop, "linfx") == 0)
      { return image_window_op_LINFX; }
    else if (strcmp(chop, "linfy") == 0)
      { return image_window_op_LINFY; }
    else if (strcmp(chop, "linvar") == 0)
      { return image_window_op_LINVAR; }
    else if (strcmp(chop, "lindev") == 0)
      { return image_window_op_LINDEV; }
    else 
      { /* Invalid op name: */
        return image_window_op_NUM;
      }
  }

const char *image_window_op_to_string(image_window_op_t op)
  {
    switch(op)
      { 
        case image_window_op_AVERAGE:    return "average";
        case image_window_op_VARIANCE:   return "variance";
        case image_window_op_DEVIATION:  return "deviation";
        case image_window_op_F:          return "f";
        case image_window_op_FX:         return "fx";
        case image_window_op_FY:         return "fy";
        case image_window_op_FXX:        return "fxx";
        case image_window_op_FXY:        return "fxy";
        case image_window_op_FYY:        return "fyy";
        case image_window_op_FXXY:       return "fxxy";
        case image_window_op_FXYY:       return "fxyy";
        case image_window_op_FXXYY:      return "fxxyy";
        case image_window_op_GRADIENT:   return "gradient";
        case image_window_op_LAPLACIAN:  return "laplacian";
        case image_window_op_ORTHICITY:  return "orthicity";
        case image_window_op_ELONGATION: return "elongation";
        case image_window_op_LINF:       return "linf";
        case image_window_op_LINFX:      return "linfx";
        case image_window_op_LINFY:      return "linfy";
        case image_window_op_LINVAR:     return "linvar";
        case image_window_op_LINDEV:     return "lindev";
        default: demand(FALSE, "invalid {op}");
      }
  }

/* SPECIFIC OPERATORS */

double image_window_op_f(int32_t ictr, int32_t nx, double smp[])
  { return smp[ictr]; }

double image_window_op_fx_unsmoothed(int32_t ictr, int32_t nx, double smp[])
  { return (smp[ictr+1] - smp[ictr-1])/2; }

double image_window_op_fx_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double sm = (smp[ictr-nx-1] + 2*smp[ictr-1] + smp[ictr+nx-1])/4;
    double sp = (smp[ictr-nx+1] + 2*smp[ictr+1] + smp[ictr+nx+1])/4;
    return (sp - sm)/2;
  }

double image_window_op_fy_unsmoothed(int32_t ictr, int32_t nx, double smp[])
  { return (smp[ictr+nx] - smp[ictr-nx])/2; }

double image_window_op_fy_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double sm = (smp[ictr-1-nx] + 2*smp[ictr-nx] + smp[ictr+1-nx])/4;
    double sp = (smp[ictr-1+nx] + 2*smp[ictr+nx] + smp[ictr+1+nx])/4;
    return (sp - sm)/2;
  }

double image_window_op_fxx_unsmoothed(int32_t ictr, int32_t nx, double smp[])
  { return smp[ictr+1] + smp[ictr-1] - 2*smp[ictr]; }

double image_window_op_fxx_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double sm = smp[ictr-nx-1] + 2*smp[ictr-1] + smp[ictr+nx-1];
    double so = smp[ictr-nx] + 2*smp[ictr] + smp[ictr+nx];
    double sp = smp[ictr-nx+1] + 2*smp[ictr+1] + smp[ictr+nx+1];
    return (sm - 2*so + sp)/4;
  }

double image_window_op_fxy(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    return (smp[ictr+nx+1] - smp[ictr+nx-1] - smp[ictr-nx+1] + smp[ictr-nx-1])/4;
  }

double image_window_op_fyy_unsmoothed(int32_t ictr, int32_t nx, double smp[])
  { return smp[ictr+nx] + smp[ictr-nx] - 2*smp[ictr]; }

double image_window_op_fyy_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double sm = smp[ictr-1-nx] + 2*smp[ictr-nx] + smp[ictr+1-nx];
    double so = smp[ictr-1] + 2*smp[ictr] + smp[ictr+1];
    double sp = smp[ictr-1+nx] + 2*smp[ictr+nx] + smp[ictr+1+nx];
    return (sm - 2*so + sp)/4;
  }

double image_window_op_fxxy(int32_t ictr, int32_t nx, double smp[])
  { double sp = smp[ictr+nx-1] - 2*smp[ictr+nx] + smp[ictr+nx+1];
    double sm = smp[ictr-nx-1] - 2*smp[ictr-nx] + smp[ictr-nx+1];
    return (sp - sm)/2;
  }

double image_window_op_fxyy(int32_t ictr, int32_t nx, double smp[])
  { double sp = smp[ictr+1-nx] - 2*smp[ictr+1] + smp[ictr+1+nx];
    double sm = smp[ictr-1-nx] - 2*smp[ictr-1] + smp[ictr-1+nx];
    return (sp - sm)/2;
  }

double image_window_op_fxxyy(int32_t ictr, int32_t nx, double smp[])
  { double sp = smp[ictr+nx-1] - 2*smp[ictr+nx] + smp[ictr+nx+1];
    double so = smp[ictr-1]    - 2*smp[ictr]    + smp[ictr+1];
    double sm = smp[ictr-nx-1] - 2*smp[ictr-nx] + smp[ictr-nx+1];
    return sp -2*so + sm;
  }

double image_window_op_gradient_unsmoothed(int32_t ictr, int32_t nx, double smp[])
  { double fx = (smp[ictr+1] - smp[ictr-1])/2;
    double fy = (smp[ictr+nx] - smp[ictr-nx])/2;
    return hypot(fx,fy);
  }

double image_window_op_gradient_unsmoothed_squared(int32_t ictr, int32_t nx, double smp[])
  { double fx = (smp[ictr+1] - smp[ictr-1])/2;
    double fy = (smp[ictr+nx] - smp[ictr-nx])/2;
    return fx*fx + fy*fy;
  }

double image_window_op_gradient_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double fx = image_window_op_fx_smoothed(ictr, nx, smp);
    double fy = image_window_op_fy_smoothed(ictr, nx, smp);
    return hypot(fx,fy);
  }

double image_window_op_gradient_smoothed_squared(int32_t ictr, int32_t nx, double smp[])
  { double fx = image_window_op_fx_smoothed(ictr, nx, smp);
    double fy = image_window_op_fy_smoothed(ictr, nx, smp);
    return fx*fx + fy*fy;
  }

double image_window_op_laplacian_unsmoothed(int32_t ictr, int32_t nx, double smp[])
  { /* Uses cross samples: */
    return smp[ictr+1] + smp[ictr-1] + smp[ictr+nx] + smp[ictr-nx] - 4*smp[ictr];
  }

double image_window_op_laplacian_smoothed(int32_t ictr, int32_t nx, double smp[])
  { /* Uses corner samples: */
    return (smp[ictr+nx+1] + smp[ictr+nx-1] + smp[ictr-nx+1] + smp[ictr-nx-1] - 4*smp[ictr])/2;
  }

double image_window_op_orthicity(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    return smp[ictr+1] + smp[ictr-1] - smp[ictr+nx] - smp[ictr-nx];
  }

double image_window_op_elongation(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    double ort = smp[ictr+1] + smp[ictr-1] - smp[ictr+nx] - smp[ictr-nx]; /* Orthicity. */
    double wrp = (smp[ictr+nx+1] - smp[ictr+nx-1] - smp[ictr-nx+1] + smp[ictr-nx-1])/2; /* Warp = {2*fxy}. */
    return hypot(wrp, ort);
  }

double image_window_op_elongation_squared(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    double ort = smp[ictr+1] + smp[ictr-1] - smp[ictr+nx] - smp[ictr-nx]; /* Orthicity. */
    double wrp = (smp[ictr+nx+1] - smp[ictr+nx-1] - smp[ictr-nx+1] + smp[ictr-nx-1])/2; /* Warp = {2*fxy}. */
    return wrp*wrp + ort*ort;
  }

double image_window_op_average(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    return 
      ( 4.0 * smp[ictr] + 
        2.0 * (smp[ictr+1] + smp[ictr-1] + smp[ictr+nx] + smp[ictr-nx]) +  
        smp[ictr+1+nx] + smp[ictr+1-nx] + smp[ictr-1+nx] + smp[ictr-1-nx]
      ) / 16.0;
  }

double image_window_op_variance(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    double avg = image_window_op_average(ictr, nx, smp);
    
    double dmm = smp[ictr-1-nx] - avg;
    double dmo = smp[ictr-1] - avg;
    double dmp = smp[ictr-1+nx] - avg;
    
    double dom = smp[ictr-nx] - avg;
    double doo = smp[ictr] - avg;
    double dop = smp[ictr+nx] - avg;
    
    double dpm = smp[ictr+1-nx] - avg;
    double dpo = smp[ictr+1] - avg;
    double dpp = smp[ictr+1+nx] - avg;
    
    return 
      ( 4.0 * doo*doo + 
        2.0 * (dmo*dmo + dpo*dpo + dom*dom + dop*dop) +  
        dmm*dmm + dmp*dmp + dpm*dpm + dpp*dpp
      ) / 16.0;
  }

double image_window_op_deviation(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    double var = image_window_op_variance(ictr, nx, smp);
    return sqrt(var);
  }

double image_window_op_linvar(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    double linf = image_window_op_average(ictr, nx, smp);
    double linfx = image_window_op_fx_smoothed(ictr, nx, smp);
    double linfy = image_window_op_fy_smoothed(ictr, nx, smp);
    
    double dmm = smp[ictr-1-nx] - (linf - linfx - linfy);
    double dmo = smp[ictr-1] - (linf - linfx);
    double dmp = smp[ictr-1+nx] - (linf - linfx + linfy) ;
    
    double dom = smp[ictr-nx] - (linf - linfy);
    double doo = smp[ictr] - linf;
    double dop = smp[ictr+nx] - (linf + linfy);
    
    double dpm = smp[ictr+1-nx] - (linf + linfx - linfy);
    double dpo = smp[ictr+1] - (linf + linfx);
    double dpp = smp[ictr+1+nx] - (linf + linfx + linfy);
    
    return 
      ( 4.0 * doo*doo + 
        2.0 * (dmo*dmo + dpo*dpo + dom*dom + dop*dop) +  
        dmm*dmm + dmp*dmp + dpm*dpm + dpp*dpp
      ) / 16.0;
  }

double image_window_op_lindev(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    double var = image_window_op_linvar(ictr, nx, smp);
    return sqrt(var);
  }

