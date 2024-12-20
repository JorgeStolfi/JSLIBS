/* See {frgb_path.h}. */
/* Last edited on 2024-12-05 07:48:33 by stolfi */ 

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <frgb.h>
#include <jsmath.h>
#include <frgb_ops.h>
#include <frgb_interp_vis.h>

#include <frgb_path.h>

#define AUTHOR \
  "  Created jan/2007 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP).\n" \
  "  The {unsigned_0} color path was developed by J. Stolfi around 2021-12" \
  " for the HotPath project with {R. Minetto}."

double frgb_path_clip_unsigned(double z);
  /* Clips {z} to the range {[0_1]}. */

double frgb_path_clip_signed(double z);
  /* Clips {z} to the range {[-1_+1]}. */

#define YR frgb_YR
#define YG frgb_YG
#define YB frgb_YB
  /* Luminance weights according to European TV standard. */

frgb_t frgb_path_map_unsigned_0(double z, int32_t cycles)
  { 
    bool_t debug = FALSE;
    
    /* Reduce the path parameter to {[0_1]}: */
    z = frgb_path_clip_unsigned(z);
    
    /* Path endpoints in HTY(UV) color space, initially for grays: */
    frgb_t q0 = (frgb_t){{ 0.0, 0.0, 0.0 }};
    frgb_t q1 = (frgb_t){{ 0.0, 0.0, 1.0 }};
    
    if (cycles != 0)
      { /* Not gray scale. Set max relative saturation {T}: */
        q0.c[1] = 1.0; q1.c[1] = 1.0;
        /* Set hue {H}: */
        double Hora = 0.2; /* Orange. */
        double Hpur = 0.8; /* Purple. */
        if (cycles > 0)
          { /* Clockwise from purple to orange: */
            q0.c[0] = (float)Hpur; q1.c[0] = (float)(Hora - abs(cycles));
          }
        else
          { /* Counterclockwise from orange to purple: */
            q0.c[0] = (float)Hora; q1.c[0] = (float)(Hpur + abs(cycles)-1);
          }
      }
    if (debug)
      { frgb_print(stderr, "q0 = ( ", &q0, 3, "%+8.4f", " )\n");
        frgb_print(stderr, "q1 = ( ", &q1, 3, "%+8.4f", " )\n");
      }
    frgb_t q = frgb_interp_vis_HTY(z, &q0, &q1, FALSE);
    frgb_from_HTY(&q);
    return q;
  }

int32_t frgb_path_map_unsigned_max_style(void)
  { return 0; }

frgb_t frgb_path_map_unsigned(double z, int32_t cycles, int32_t style)
  { switch(style)
      { case 0: return frgb_path_map_unsigned_0(z, cycles);
        default: demand(FALSE, "invalid {style}");
      }
  }

frgb_t frgb_path_map_signed_0(double z, int32_t cycles)
  { /* A simple S-shaped path in the cube.  Ignores the {cycles} argument. */
    /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Choose extrmal colors: */
    float cmin[3] = { 0.500f, 0.500f, 0.500f }; /* Middle gray. */
    float cmed[3] = { 1.000f, 0.333f, 0.000f }; /* B�zier middle color. */
    float cmax[3] = { 1.000f, 0.933f, 0.900f }; /* Color for max positive spots. */
    /* Interpolate between them with {z}: */
    frgb_t clr;
    if (fabs(z) <= 1.0e-5)
      { /* Map to center gray: */
        for (int32_t kc = 0;  kc < 3; kc++) { clr.c[kc] = 0.500; }
      }
    else
      { /* Stretch argument away from 0: */
        double t = sqrt(fabs(z));
        /* Interpolate between {cmin} and {cmax} with {t} as ratio: */
        double s = 1 - t;
        for (int32_t kc = 0;  kc < 3; kc++) 
          { /* DeCasteljau's quadratic interpolation algorithm: */
            float v0t = (float)(s*cmin[kc] + t*cmed[kc]);
            float vt1 = (float)(s*cmed[kc] + t*cmax[kc]);
            float vtt = (float)(s*v0t + t*vt1);
            /* If {z} is negative, complement with respect to center  gray: */
            if (z < 0) { vtt = 1.000f - vtt; }
            clr.c[kc] = vtt;
          }
      }
    return clr;
  }

frgb_t frgb_path_map_signed_1(double z, int32_t cycles)
  { /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Compute {z} point along a spiral path: */
    frgb_t clr;
    if (fabs(z) == 0.0)
      { /* Map to center gray: */
        for (int32_t kc = 0;  kc < 3; kc++) { clr.c[kc] = 0.500; }
      }
    else
      { /* Luminance of min and max colors: */
        /* double ymax = + YR*cmax[0] + YG*cmax[1] + YB*cmax[2]; */
        /* double ymin = + YR*cmin[0] + YG*cmin[1] + YB*cmin[2]; */
        double ymax = 0.8667;
        double ymin = 0.5333;

        /* The luminance interpolates from {ymax} to {ymin} with {abs(z)} as ratio: */
        double ah = fabs(z);
        double smax = ah, smin = 1.0 - ah;
        double y = smin*ymin + smax*ymax;    
      
        /* Middle gray, and two color displacements of zero luminance: */
        float m[3] = { +0.5000f, +0.5000f, +0.5000f };
        float u[3] = { +0.5000f, -0.2500f, 00.0000f };
        float v[3] = { 00.0000f, +0.0833f, -0.5000f };
        /* Note that {m + au*u + av*v} is valid for any {au,av} in {[-1_+1]}. */

        /* Map {abs(z)} to a point {au,av} on the {[-1_+1]�[0_+1]} square: */
        /* Use an ellipse inscribed in that rectangle, with {cycles} turns: */
        double satmin = 0.1; /* Min saturation of colors on ellipse. */
        double phase = ah*cycles*2*M_PI; /* Phase on ellipse. */
        double uctr = (1 + satmin)/2, vctr = uctr; /* Center or ellipse. */
        double urad = (1 - satmin)/2, vrad = urad; /* Semidiameters of ellipse. */
        double au = uctr + urad*cos(phase);
        double av = vctr + vrad*sin(phase);

        /* Make {A} a mix of {m + au*u + av*v} with white, ratio {y-0.5}: */
        double qwht = y - 0.5000, qell = 1 - qwht;
        /* A = qwht*(1,1,1) + qell*(m + au*u + av*v); */

        /* Make {B} a mix of {m} and {u+v}, ratio {ah}: */
        double puv = ah, pm = 1 - ah;
        /* B = pm*m + puv*u + puv*v; */

        /* Make a a mix {A} and {B}, with {B} significant only for small {ah}: */
        double eps = 0.1; /* Mixing starts effectively below this. */
        double k = pow(eps/ah, 2);
        double klin = k/(1 + k), kspr = 1/(1 + k);
        /* C = kspr*A + klin*B; */

        for (int32_t kc = 0;  kc < 3; kc++) 
          { double Ac = qwht*1.0 + qell*(m[kc] + au*u[kc] + av*v[kc]);
            double Bc = pm*m[kc] + puv*u[kc] + puv*v[kc];
            double Cc = kspr*Ac + klin*Bc;
            /* If {z} is negative, complement with respect to center  gray: */
            if (z < 0) { Cc = 1.000 - Cc; }
            clr.c[kc] = (float)Cc;
          }
      }
    
    return clr;
  }

frgb_t frgb_path_map_signed_2(double z, int32_t cycles)
  { /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Compute {z} point along a spiral path: */
    frgb_t clr;
    if (fabs(z) == 0.0)
      { /* Map to center gray: */
        for (int32_t kc = 0;  kc < 3; kc++) { clr.c[kc] = 0.500; }
      }
    else
      { /* Luminance of min and max colors: */
        /* double ymax = + YR*cmax[0] + YG*cmax[1] + YB*cmax[2]; */
        /* double ymin = + YR*cmin[0] + YG*cmin[1] + YB*cmin[2]; */
        double ymax = 0.8667;
        double ymin = 0.5333;

        /* The luminance interpolates from {ymax} to {ymin} with {abs(z)} as ratio: */
        double ah = fabs(z);
        double smax = ah, smin = 1.0 - ah;
        double y = smin*ymin + smax*ymax;    
      
        /* Middle gray, and two color displacements of zero luminance: */
        float m[3] = { +0.5000f, +0.5000f, +0.5000f };
        float u[3] = { +0.5000f, -0.2500f, 00.0000f };
        float v[3] = { 00.0000f, +0.0833f, -0.5000f };
        /* Note that {m + au*u + av*v} is valid for any {au,av} in {[-1_+1]}. */

        /* Map {abs(z)} to a point {au,av} on the {[-1_+1]�[-1_+1]} square: */
        /* Use an ellipse inscribed in that rectangle, with {cycles} turns: */
        double phase = ah*cycles*2*M_PI; /* Phase on ellipse. */
        double uctr = 0.0, vctr = 0.0; /* Center or ellipse. */
        double urad = 1.0, vrad = 1.0; /* Semidiameters of ellipse. */
        double au = uctr + urad*cos(phase);
        double av = vctr + vrad*sin(phase);

        /* Make {A} a mix of {m + au*u + av*v} with white, ratio {y-0.5}: */
        double qwht = y - 0.5000, qell = 1 - qwht;
        /* A = qwht*(1,1,1) + qell*(m + au*u + av*v); */

        /* Make a a mix {A} and {m}, with {m} significant only for small {ah}: */
        double eps = 0.1; /* Mixing starts effectively below this. */
        double k = pow(eps/ah, 2);
        double kmid = k/(1 + k), kspr = 1/(1 + k);
        /* C = kspr*A + kmid*B; */

        for (int32_t kc = 0;  kc < 3; kc++) 
          { double Ac = qwht*1.0 + qell*(m[kc] + au*u[kc] + av*v[kc]);
            double Cc = kspr*Ac + kmid*m[kc];
            /* If {z} is negative, complement with respect to center  gray: */
            if (z < 0) { Cc = 1.000 - Cc; }
            clr.c[kc] = (float)Cc;
          }
      }
    
    return clr;
  }
    
int32_t frgb_path_map_signed_max_style(void)
  { return 2; }
  
frgb_t frgb_path_map_signed(double z, int32_t cycles, int32_t style)
  { 
    switch(style) 
      {
        case 0: return frgb_path_map_signed_0(z, cycles); break;
        case 1: return frgb_path_map_signed_1(z, cycles); break;
        case 2: return frgb_path_map_signed_2(z, cycles); break;
        default: demand(FALSE, "invalid {style}"); 
      }
  }

double frgb_path_reduce_unsigned(double z)
  { if ((z < 0) || (z > 1))
      { /* Compute {fz = z} reduced modulo 1: */
        double fz = z - floor(z);
        /* Make sure {fz} is in {[0_1]}: */
        while (fz < 0.0) { fz += 1.0; }
        while (fz > 1.0) { fz -= 1.0; }
        /* If {z} was positive, keep {fz} away from 0: */
        if ((z > 0) && (fz == 0)) { fz = 1.0; }
        /* If {z} was negative, keep {fz} away from 1: */
        if ((z <= 0) && (fz == 1)) { fz = 1.0; }
        return fz;
      }
    else
      { return z; }
  }

double frgb_path_reduce_signed(double z)
  { if ((z < -1) || (z > +1))
      { /* Compute the signed fractional part {fz} of {z}: */
        double q = floor(fabs(z));
        double fz = z - (z < 0 ? -q : +q);
        /* Make sure {fz} is in {[-1_+1]}: */
        while (fz < -1.0) { fz += 1.0; }
        while (fz > +1.0) { fz -= 1.0; }
        /* If {z} was nonzero, make {fz} so: */
        if ((z < 0) && (fz == 0)) { fz = -1.0; }
        if ((z > 0) && (fz == 0)) { fz = +1.0; }
        return fz;
      }
    else
      { return z; }
  }
    
double frgb_path_clip_unsigned(double z)
  { if (z < 0) 
      { return 0; }
    else if (z > 1)
      { return 1; }
    else
      { return z; }
  }

double frgb_path_clip_signed(double z)
  { if (z < -1) 
      { return -1; }
    else if (z > +1)
      { return +1; }
    else
      { return z; }
  }

