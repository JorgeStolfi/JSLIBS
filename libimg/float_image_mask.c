/* See {float_image_mask.h}. */
/* Last edited on 2013-10-21 00:05:41 by stolfilocal */

#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
 
#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_mask.h>

double float_image_mask_profile(double iu2, double ou2, int ord);
  /* Evaluates a uni-dimensional windowing function {F(u2)}
    with inner parameter {iu2}, outer paramter {ou2}, and
    order {ord}, at argument {u2 = 1.0}. */

void float_image_mask_window(float_image_t *msk, int ic, int ord, bool_t round)
  { /* Get the window dims: */
    int NC = (int)msk->sz[1]; if (NC == 0) { return; }
    int NX = (int)msk->sz[1]; if (NX == 0) { return; }
    int NY = (int)msk->sz[2]; if (NY == 0) { return; }
    
    /* Argument checking: */
    demand((ic >= 0) && (ic < NC), "invalid channel");
    demand((ord >= 0) && (ord <= 2), "invalid {ord}");
    
    /* Compute the domain's center {xc,yc}: */
    double xc = ((double)NX)/2;
    double yc = ((double)NY)/2;

    /* Compute inner semidiameters {irx,iry}. */
    double margin = 2 - ord;
    double irx = fmax(xc / 2, xc - margin);
    double iry = fmax(yc / 2, yc - margin);
    
    /* Compute outer semidiameters {orx,ory}: */
    double orx = xc;
    double ory = yc;
    
    /* Enumerate all pixels. */
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { /* Pixel center coords relative to window center: */
            double xp = ix + 0.5 - xc;
            double yp = iy + 0.5 - yc;
            
            /* Ratios of {xp,yp} to inner semidiameters, squared: */
            double iux2 = xp/irx; iux2 *= iux2;
            double iuy2 = yp/iry; iuy2 *= iuy2;
            
            /* Ratios of {xp,yp} to outer semidiameters, squared: */
            double oux2 = xp/orx; oux2 *= oux2; 
            double ouy2 = yp/ory; ouy2 *= ouy2;
            
            /* Compute antialiased value {val}: */
            double val;
            if (round)
              { val = float_image_mask_profile(iux2+iuy2, oux2+ouy2, ord); }
            else
              { double valx =  float_image_mask_profile(iux2, oux2, ord);
                double valy =  float_image_mask_profile(iuy2, ouy2, ord);
                val = valx * valy;
              }
              
            float_image_fill_pixel(msk, ix, iy, (float)val);
          }
      }
  }
  
double float_image_mask_profile(double iu2, double ou2, int ord)
  {
    demand(iu2 >= ou2, "bad {iu2,ou2}");
    if (ou2 >= 1.0) { /* Fully outside: */ return 0; }
    switch (ord)
      {
      case 0: 
        { if (iu2 <= 1.0) 
            { /* Fully inside: */ 
              return 1.0;
            } 
          else
            { /* Smooth cubic sigmoid: */ 
              double d = iu2 - ou2;
              double s = (1 - ou2)/d;
              return s*s*(3 - 2*s);
            }
        }
      case 1: 
        { double mu2 = (iu2 + ou2)/2;
          if (iu2 <= 1.0) 
            { /* Fully inside: */ 
              return 1.0 - mu2;
            } 
          else
            { /* Smooth parabolic fillet: */ 
              double d = iu2 - ou2;
              double a = 1.0 - ou2;
              return a*a/d/2;
            }
        }
      case 2:
        { demand(iu2 == ou2, "bad {iu2,ou2} for {ord=2}");
          /* Biquadratic Hann: */
          double d = 1 - ou2;
          return d*d;
        }
      default:
        demand(FALSE, "invalid {ord}");
      }
  }

void float_image_mask_mul_gauss
  ( float_image_t *msk, 
    int ic, 
    double sx, 
    double sy
  )
  { /* Get the window dims: */
    int NC = (int)msk->sz[1];
    int NX = (int)msk->sz[1];
    int NY = (int)msk->sz[2];
    
    /* Argument checking: */
    demand((ic >= 0) && (ic < NC), "invalid channel");
    demand(sx > 0, "bad {sx}");
    demand(sy > 0, "bad {sy}");
    
    /* Compute the domain's center {xc,yc}: */
    double xc = ((double)NX)/2;
    double yc = ((double)NY)/2;
     
    /* Multiply all pixels by {F(x,y)}. */
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { float *p = float_image_get_sample_address(msk, ic, ix, iy);
            double wxy = (*p);
            if (wxy != 0.0)
              { double tx = (ix + 0.5 - xc)/sx;
                double ty = (iy + 0.5 - yc)/sy;
                double r2 = tx*tx + ty*ty;
                double F = exp(-r2/2);
                (*p) = (float)(F * wxy);
              }
          }
      }
  }

void float_image_mask_mul_power
  ( float_image_t *msk, 
    int ic, 
    double sx, 
    double sy, 
    double pwr
  )
  { /* Get the window dims: */
    int NC = (int)msk->sz[1];
    int NX = (int)msk->sz[1];
    int NY = (int)msk->sz[2];
    
    /* Argument checking: */
    demand((ic >= 0) && (ic < NC), "invalid channel");
    demand(sx > 0, "bad {sx}");
    demand(sy > 0, "bad {sy}");
    
    /* Compute the domain's center {xc,yc}: */
    double xc = ((double)NX)/2;
    double yc = ((double)NY)/2;
     
    double e = -0.5*pwr;
    
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { float *p = float_image_get_sample_address(msk, ic, ix, iy);
            double wxy = (*p);
            if (wxy != 0.0)
              { double tx = (ix + 0.5 - xc)/sx;
                double ty = (iy + 0.5 - yc)/sy;
                double r2 = tx*tx + ty*ty + 1.0;
                double F = pow(r2, e);
                (*p) = (float)(F * wxy);
              }
          }
      }
  }


float_image_mask_stats_t float_image_mask_stats_get(float_image_t *msk, int ic)
  {
    int NC = (int)msk->sz[0];
    int NX = (int)msk->sz[1];
    int NY = (int)msk->sz[2];
    
    demand((ic >= 0) && (ic < NC), "invalid channel"); 

    /* Compute the domain's center {xc,yc}: */
    double xc = ((double)NX)/2;
    double yc = ((double)NY)/2;

    /* Return result: */
    float_image_mask_stats_t S;
    
    /* Initialize counts and ranges: */
    S.NX = NX;
    S.NY = NY;
    S.nNAN = 0;
    S.nINF = 0;
    S.min = +INF;
    S.max = -INF;
    S.ext[0] = S.ext[1] = 0.0;
    S.rad = 0;
    
    int nOK = 0;            /* Count of samples that are neither NAN nor {±INF}. */
    double sum_v = 0;       /* Sum of image values {v}. */
    double sum_v2 = 0;      /* Sum of squared image values {v^2}. */
    double sum_av = 0;      /* Sum of {|v|}. */
    double sum_av_x = 0;    /* Sum of {|v|*x}. */
    double sum_av_y = 0;    /* Sum of {|v|*y}. */
    double sum_av_dx2 = 0;  /* Sum of {|v|*(x - xc)^2}. */
    double sum_av_dy2 = 0;  /* Sum of {|v|*(y - yc)^2}. */
    
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { double v = float_image_get_sample(msk, ic, ix, iy);
            if (isnan(v))
              { /* Tally {NAN}s: */  
                S.nNAN++;
              }
            else if (isinf(v) != 0)
              { /* Tally {±INF}s: */
                S.nINF++;
              }
            else
              { /* Tally finite samples: */
                nOK++;
                /* Uptate the sample range: */
                S.min = fmin(S.min, v);
                S.max = fmax(S.max, v);
                /* Accumulate sample sums: */
                sum_v += v;
                sum_v2 += v*v;
                /* Accumulate sample barycenter and moment sums: */
                double av = fabs(v);
                sum_av += av;
                double x = ix + 0.5;
                double y = iy + 0.5;
                sum_av_x += av*x;
                sum_av_y += av*y;
                double dx = fabs(x - xc);
                double dy = fabs(y - yc);
                /* Assume that weight is uniformly distributed in pixel: */
                sum_av_dx2 += av*(dx*dx + 1.0/12.0);
                sum_av_dy2 += av*(dy*dy + 1.0/12.0);
                if (v != 0)
                  { /* Update sample extents and radius: */
                    S.ext[0] = fmax(S.ext[0], dx + 0.5);
                    S.ext[1] = fmax(S.ext[1], dy + 0.5);
                    S.rad = fmax(S.rad, hypot(dx+0.5, dy+0.5));
                  }
              }
          }
      }
    
    /* Compute sample averages and deviations: */
    S.avg = sum_v/nOK;
    S.dev = sqrt(fmax(0.0, sum_v2/nOK - S.avg*S.avg));
    
    /* Compute sample barycenter: */
    S.ctr[0] = sum_av_x/sum_av;
    S.ctr[1] = sum_av_y/sum_av;
            
    /* Compute sample moments: */
    S.mmt[0] = sqrt(sum_av_dx2/sum_av);
    S.mmt[1] = sqrt(sum_av_dy2/sum_av);
    
    return S;
  }

void float_image_mask_stats_print(FILE *wr, float_image_mask_stats_t *S)
  {
    int NX = S->NX;
    int NY = S->NY;
    
    /* Compute the domain's center {xc,yc}: */
    double xc = ((double)NX)/2;
    double yc = ((double)NY)/2;

    fprintf(wr, "  dimensions =        ( %5d %5d )\n", S->NX, S->NY);
    fprintf(wr, "  total samples =     %9d\n", S->NX*S->NY);
    fprintf(wr, "  infinite samples =  %9d\n", S->nINF);
    fprintf(wr, "  undefined samples = %9d\n", S->nNAN);
    fprintf(wr, "  sample average =    %10.4f\n", S->avg);
    fprintf(wr, "  sample deviation =  %10.4f\n", S->dev);
    fprintf(wr, "  sample range =      [ %8.4f %8.4f ]\n", S->min, S->max);
    fprintf(wr, "  geometric center =  ( %8.3f %8.3f )\n", xc, yc);
    fprintf(wr, "  geometric extents = ( %8.3f %8.3f )\n", xc, yc);
    fprintf(wr, "  geometric radius =  %10.3f\n", hypot(xc,yc));
    fprintf(wr, "  barycenter =        ( %8.3f %8.3f )\n", S->ctr[0], S->ctr[1]);
    fprintf(wr, "  actual extents =    ( %8.3f %8.3f )\n", S->ext[0], S->ext[1]);
    fprintf(wr, "  actual radius =     %10.3f\n", S->rad);
    fprintf(wr, "  mass moments =      ( %8.3f %8.3f )\n", S->mmt[0], S->mmt[1]);
    fprintf(wr, "  total moment =      %10.3f\n", hypot(S->mmt[0], S->mmt[1]));
  }
