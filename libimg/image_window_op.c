/* See {image_window_op.h}. */
/* Last edited on 2020-11-15 16:16:36 by jstolfi */

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
            case image_window_op_IDENT:      nx = 1; ny = 1; break;
            case image_window_op_DX:         nx = 3; ny = 1; break;
            case image_window_op_DY:         nx = 1; ny = 3; break;
            case image_window_op_DXX:        nx = 3; ny = 1; break;
            case image_window_op_DXY:        nx = 3; ny = 3; break;
            case image_window_op_DYY:        nx = 1; ny = 3; break;
            case image_window_op_GRADIENT:   nx = 3; ny = 3; break;
            case image_window_op_LAPLACIAN:  nx = 3; ny = 3; break;
            case image_window_op_ORTHICITY:  nx = 3; ny = 3; break;
            case image_window_op_ELONGATION: nx = 3; ny = 3; break;
            case image_window_op_AVERAGE:    nx = 3; ny = 3; break;
            case image_window_op_DEVIATION:  nx = 3; ny = 3; break;
            default: demand(FALSE, "invalid {op}"); 
          }
      }
    (*nxP) = nx;
    (*nyP) = ny;
    assert(nx % 2 == 1);
    assert(ny % 2 == 1);
    (*ictrP) = (ny/2)*nx + (nx/2);
  }

void image_window_op_get_range(image_window_op_t op, bool_t squared, double *loP, double *hiP)
  { 
    /* Determine the natural range {[lo_hi]} of the operator, without squaring: */
    double lo = NAN, hi = NAN;
    if (! squared)
      {
        switch(op)
          {
            case image_window_op_IDENT:      lo = 00.0; hi = +1.0;       break;
            case image_window_op_DX:         lo = -0.5; hi = +0.5;       break;
            case image_window_op_DY:         lo = -0.5; hi = +0.5;       break;
            case image_window_op_DXX:        lo = -2.0; hi = +2.0;       break;
            case image_window_op_DXY:        lo = -0.5; hi = +0.5;       break;
            case image_window_op_DYY:        lo = -0.5; hi = +0.5;       break;
            case image_window_op_GRADIENT:   lo = 00.0; hi = +sqrt(0.5); break;
            case image_window_op_LAPLACIAN:  lo = -4.0; hi = +4.0;       break;
            case image_window_op_ORTHICITY:  lo = -2.0; hi = +2.0;       break;
            case image_window_op_ELONGATION: lo = 00.0; hi = +sqrt(5.0); break;
            case image_window_op_AVERAGE:    lo = 00.0; hi = +1.0;       break;
            case image_window_op_DEVIATION:  lo = 00.0; hi = +0.5;       break;
            default: demand(FALSE, "invalid {op}"); 
          }
      }
    else
      { lo = 0.0;
        switch(op)
          {
            case image_window_op_IDENT:      hi = +1.00; break;
            case image_window_op_DX:         hi = +0.25; break;
            case image_window_op_DY:         hi = +0.25; break;
            case image_window_op_DXX:        hi = +4.00; break;
            case image_window_op_DXY:        hi = +0.25; break;
            case image_window_op_DYY:        hi = +0.50; break;
            case image_window_op_GRADIENT:   hi = +0.50; break;
            case image_window_op_LAPLACIAN:  hi = +16.0; break;
            case image_window_op_ORTHICITY:  hi = +4.00; break;
            case image_window_op_ELONGATION: hi = +5.00; break;
            case image_window_op_AVERAGE:    hi = +1.00; break;
            case image_window_op_DEVIATION:  hi = +0.25; break;
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
      { case image_window_op_IDENT:
          if (smoothed)
            { res = image_window_op_average(ictr, nx, smp); }
          else
            { res = smp[ictr]; }
          if (squared) { res = res*res; }
          break;
        case image_window_op_DX:
          if (smoothed)
            { res = image_window_op_dx_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_dx(ictr, nx, smp); }
          if (squared) { res = res*res; }
          break;
        case image_window_op_DY:
          if (smoothed)
            { res = image_window_op_dy_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_dy(ictr, nx, smp); }
          if (squared) { res = res*res; }
          break;
        case image_window_op_DXX:
          if (smoothed)
            { res = image_window_op_dxx_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_dxx(ictr, nx, smp); }
          if (squared) { res = res*res; }
          break;
        case image_window_op_DXY:
          res = image_window_op_dxy(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_DYY:
          if (smoothed)
            { res = image_window_op_dyy_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_dyy(ictr, nx, smp); }
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
                { res = image_window_op_gradient_squared(ictr, nx, smp); }
              else
                { res = image_window_op_gradient(ictr, nx, smp); }
            }
          break;
        case image_window_op_LAPLACIAN:
          if (smoothed)
            { res = image_window_op_laplacian_smoothed(ictr, nx, smp); }
          else
            { res = image_window_op_laplacian(ictr, nx, smp); }
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
        case image_window_op_AVERAGE:
          res = image_window_op_average(ictr, nx, smp);
          if (squared) { res = res*res; }
          break;
        case image_window_op_DEVIATION:
          if (squared) 
            { res = image_window_op_deviation_squared(ictr, nx, smp); }
          else
            { res = image_window_op_deviation(ictr, nx, smp); }
          break;
        default: demand(FALSE, "invalid {op}"); 
      }
    return res;
  }

image_window_op_t image_window_op_from_string(const char *chop)
  {
    if (strcmp(chop, "ident") == 0)
      { return image_window_op_IDENT; }
    else if (strcmp(chop, "dx") == 0)
      { return image_window_op_DX; }
    else if (strcmp(chop, "dy") == 0)
      { return image_window_op_DY; }
    else if (strcmp(chop, "dxx") == 0)
      { return image_window_op_DXX; }
    else if (strcmp(chop, "dxy") == 0)
      { return image_window_op_DXY; }
    else if (strcmp(chop, "dyy") == 0)
      { return image_window_op_DYY; }
    else if (strcmp(chop, "gradient") == 0)
      { return image_window_op_GRADIENT; }
    else if (strcmp(chop, "laplacian") == 0)
      { return image_window_op_LAPLACIAN; }
    else if (strcmp(chop, "orthicity") == 0)
      { return image_window_op_ORTHICITY; }
    else if (strcmp(chop, "elongation") == 0)
      { return image_window_op_ELONGATION; }
    else if (strcmp(chop, "average") == 0)
      { return image_window_op_AVERAGE; }
    else if (strcmp(chop, "deviation") == 0)
      { return image_window_op_DEVIATION; }
    else 
      { /* Invalid op name: */
        return image_window_op_NUM_VALUES;
      }
  }

const char *image_window_op_to_string(image_window_op_t op)
  {
    switch(op)
      { 
        case image_window_op_IDENT:      return "ident";
        case image_window_op_DX:         return "dx";
        case image_window_op_DY:         return "dy";
        case image_window_op_DXX:        return "dxx";
        case image_window_op_DXY:        return "dxy";
        case image_window_op_DYY:        return "dyy";
        case image_window_op_GRADIENT:   return "gradient";
        case image_window_op_LAPLACIAN:  return "laplacian";
        case image_window_op_ORTHICITY:  return "orthicity";
        case image_window_op_ELONGATION: return "elongation";
        case image_window_op_AVERAGE:    return "average";
        case image_window_op_DEVIATION:  return "deviation";
        default: demand(FALSE, "invalid {op}");
      }
  }

/* SPECIFIC OPERATORS */

double image_window_op_ident(int32_t ictr, int32_t nx, double smp[])
  { return smp[ictr]; }

double image_window_op_dx(int32_t ictr, int32_t nx, double smp[])
  { return (smp[ictr+1] - smp[ictr-1])/2; }

double image_window_op_dx_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double sm = (smp[ictr-nx-1] + 2*smp[ictr-1] + smp[ictr+nx-1])/4;
    double sp = (smp[ictr-nx+1] + 2*smp[ictr+1] + smp[ictr+nx+1])/4;
    return (sp - sm)/2;
  }

double image_window_op_dy(int32_t ictr, int32_t nx, double smp[])
  { return (smp[ictr+nx] - smp[ictr-nx])/2; }

double image_window_op_dy_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double sm = (smp[ictr-1-nx] + 2*smp[ictr-nx] + smp[ictr+1-nx])/4;
    double sp = (smp[ictr-1+nx] + 2*smp[ictr+nx] + smp[ictr+1+nx])/4;
    return (sp - sm)/2;
  }

double image_window_op_dxx(int32_t ictr, int32_t nx, double smp[])
  { return smp[ictr+1] + smp[ictr-1] - 2*smp[ictr]; }

double image_window_op_dxx_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double sm = (smp[ictr-nx-1] + 2*smp[ictr-1] + smp[ictr+nx-1])/4;
    double so = (smp[ictr-nx] + 2*smp[ictr] + smp[ictr+nx])/4;
    double sp = (smp[ictr-nx+1] + 2*smp[ictr+1] + smp[ictr+nx+1])/4;
    return sm + sp - 2*so;
  }

double image_window_op_dxy(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    return (smp[ictr+nx+1] - smp[ictr+nx-1] - smp[ictr-nx+1] + smp[ictr-nx-1])/4;
  }

double image_window_op_dyy(int32_t ictr, int32_t nx, double smp[])
  { return smp[ictr+nx] + smp[ictr-nx] - 2*smp[ictr]; }

double image_window_op_dyy_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double sm = (smp[ictr-1-nx] + 2*smp[ictr-nx] + smp[ictr+1-nx])/4;
    double so = (smp[ictr-1] + 2*smp[ictr] + smp[ictr+1])/4;
    double sp = (smp[ictr-1+nx] + 2*smp[ictr+nx] + smp[ictr+1+nx])/4;
    return sm + sp - 2*so;
  }

double image_window_op_gradient(int32_t ictr, int32_t nx, double smp[])
  { double dx = (smp[ictr+1] - smp[ictr-1])/2;
    double dy = (smp[ictr+nx] - smp[ictr-nx])/2;
    return hypot(dx,dy);
  }

double image_window_op_gradient_squared(int32_t ictr, int32_t nx, double smp[])
  { double dx = (smp[ictr+1] - smp[ictr-1])/2;
    double dy = (smp[ictr+nx] - smp[ictr-nx])/2;
    return dx*dx + dy*dy;
  }

double image_window_op_gradient_smoothed(int32_t ictr, int32_t nx, double smp[])
  { double dx = image_window_op_dx_smoothed(ictr, nx, smp);
    double dy = image_window_op_dy_smoothed(ictr, nx, smp);
    return hypot(dx,dy);
  }

double image_window_op_gradient_smoothed_squared(int32_t ictr, int32_t nx, double smp[])
  { double dx = image_window_op_dx_smoothed(ictr, nx, smp);
    double dy = image_window_op_dy_smoothed(ictr, nx, smp);
    return dx*dx + dy*dy;
  }

double image_window_op_laplacian(int32_t ictr, int32_t nx, double smp[])
  { return smp[ictr+1] + smp[ictr-1] + smp[ictr+nx] + smp[ictr-nx] - 4*smp[ictr]; }

double image_window_op_laplacian_smoothed(int32_t ictr, int32_t nx, double smp[])
  { return (smp[ictr+1+nx] + smp[ictr-1+nx] + smp[ictr+1-nx] + smp[ictr-1-nx] - 4*smp[ictr])/2; }

double image_window_op_orthicity(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    return smp[ictr+1] + smp[ictr-1] - smp[ictr+nx] - smp[ictr-nx];
  }

double image_window_op_elongation(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    double ort = smp[ictr+1] + smp[ictr-1] - smp[ictr+nx] - smp[ictr-nx]; /* Orthicity. */
    double wrp = (smp[ictr+nx+1] - smp[ictr+nx-1] - smp[ictr-nx+1] + smp[ictr-nx-1])/2; /* Warp = {2*dxy}. */
    return hypot(wrp, ort);
  }

double image_window_op_elongation_squared(int32_t ictr, int32_t nx, double smp[])
  { /* Smoothed version is the same as unsmoothed one: */
    double ort = smp[ictr+1] + smp[ictr-1] - smp[ictr+nx] - smp[ictr-nx]; /* Orthicity. */
    double wrp = (smp[ictr+nx+1] - smp[ictr+nx-1] - smp[ictr-nx+1] + smp[ictr-nx-1])/2; /* Warp = {2*dxy}. */
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

double image_window_op_deviation_squared(int32_t ictr, int32_t nx, double smp[])
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
    double var = image_window_op_deviation_squared(ictr, nx, smp);
    return sqrt(var);
  }

