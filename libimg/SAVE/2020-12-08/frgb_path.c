/* See {frgb_path.h}. */
/* Last edited on 2023-03-06 19:54:56 by stolfi */ 

#define _GNU_SOURCE
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <frgb_path.h>

#define AUTHOR \
  "  Created jan/2007 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP).\n" \
  "  The {unsigned_0} color path function was originally developed" \
  " and implemented by J. Stolfi in nov/1997 for {recluster}, a" \
  " vector clustering program developed for Voynich manuscript" \
  " analysis.\n" \
  "  The {signed_0}, {signed_1}, and {signed_2} color path functions" \
  " were originally written by J. Stolfi in apr/2006 for {fni_view}, a" \
  " terrain viewer developed as part of the UFF-UNICAMP Photometric" \
  " Stereo project.\n" \
  "  The {signed_U} color path function was developed by J. Stolfi" \
  " in 08/jan/2007 for use in {spr2D}, a Maxwell equations integrator" \
  " being developed by UNICAMP-INPE-UPORTO, and implemented" \
  " as the {gawk} script {make-symmetric-palette}.  The variants" \
  " {unsigned_S} and {signed_S} were deeloped a few days later."

double frgb_path_clip_unsigned(double z);
  /* Clips {z} to the range {[0_1]}. */

double frgb_path_clip_signed(double z);
  /* Clips {z} to the range {[-1_+1]}. */

#define YR  (0.298911)
#define YG  (0.586611)
#define YB  (0.114478)
  /* Luminance weights according to European TV standard. */

frgb_t frgb_path_map_unsigned_0(double z)
  { 
    /* This path is based on an ad-hoc formula that produces a smooth
      (C infinity)sequence of colors, with fairly varied hues and
      monotonically increasing {Y}. Unfortunately the formula is not
      easy to change, so the {cycle} parameter is ignored.
      
      The {R} and {B} channels are computed by sinusoidal functions
      {(1-cos(k*PI*z))/2} with frequencies {k} 3 and 7, respectively.
      
      The combination of the two odd-frequency sinusoids is a
      Lissajous-like curve on the {R-B} plane, that starts at the
      {(0,0)} corner of the unit square and ends at the {(1,1)}
      corner. This curve touches the two sides {R=0} and {R=1}, four
      times each, (including the ends) and the {B=0} and {B=1} sides,
      twice each.
      
      The {G} channels is then computed so that the brightness {Y} is
      equal to {z}: namely, {G = (z - YR*R - YB*B)/YG} where
      {YR,YG,YB} are the {Y} values of pure red, pure green, and pure
      blue.
       
      By itself, this formula will yield {G<0} at the beginning of the
      curve, and {G>1} at the end. To fix these problems, the argument
      {z} in the sinusoid that defines {R} is replaced by a non-linear
      (but monotonic) map from [0_1] onto [0_1], namely {z^2*(3-2*z)}.
      This change (which does not affect the qualitative character of
      the {R-B} curve) reduces the {R} component near the beginning of
      the curve, and supresses the negative {G}s almost entirely. This
      change also fixes the symmetric {G>1} problem at the other end.
      To conclude the therapy, the {R} component is mixed with a small
      amount (3%) of {z}. (The path then fails to touch the {R=0} and
      {R=1} planes, but still gets quite close to them.) With these
      changes to the formula, {G} remains within [0_1]. Thus no
      cube-clipping is necessary, and the complete {R,G,B} path
      remains C-infinity.
      
      As a byproduct, the non-linear remapping of the {R} function
      also makes the path nearly tangent to the {R=0} plane at the
      beginning, and to the {R=1} plane at the end. */

    /* Reduce the path parameter to {[0_1]}: */
    z = frgb_path_clip_unsigned(z);

    /* Choose the red and blue half-period counts {nR,nB}: */
    int nR = 3, nB = 7;  /* Magic values... */

    /* Compute the red and blue clocks {tR,tB}: */
    double tR = z*z*(3-2*z), tB = z;

    /* Compute {R} and {B} as a Lissajous-like path in {[0_1]^2}: */
    double R = (1 - cos(nR*M_PI*tR))/2;
    double B = (1 - cos(nB*M_PI*tB))/2;
    
    /* Adjust {R} so that the next command yields {G} in {[0_1]}: */
    R = 0.03*z + 0.97*R;
    
    /* Compute {G} so that the brightness (European TV {Y}) is {z}: */
    double G = (z - YR*R - YB*B)/YG;
    
    assert((R >= 0) && (R <= 1));
    assert((G >= 0) && (G <= 1));
    assert((B >= 0) && (B <= 1));
    
    return (frgb_t){{ (float)R, (float)G, (float)B }};
  }

frgb_t frgb_path_map_unsigned_1(double z, double H0, double Y0, double H1, double Y1)
  { 
    /* Reduce the path parameter to {[0_1]}: */
    z = frgb_path_clip_unsigned(z);
    
    double tH = (1-z)*H0 + z*H1; /* Desired hue. */
    double tT = 1.0;             /* Desired relative saturation. */
    double tY = (1-z)*Y0 + z*Y1; /* Desired luminosity. */
    
    /* Compute the color {v} with hue {tH}, lum {tY} and relsat {tS}: */
    frgb_t clr = (frgb_t){{ (float)tH, (float)tT, (float)tY }};
    frgb_from_HTY_UV(&clr);
    return clr;
  }

frgb_t frgb_path_map_signed_0(double z, int cycles)
  { /* A simple S-shaped path in the cube.  Ignores the {cycles} argument. */
    /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Choose extrmal colors: */
    float cmin[3] = { 0.500f, 0.500f, 0.500f }; /* Middle gray. */
    float cmed[3] = { 1.000f, 0.333f, 0.000f }; /* B�zier middle color. */
    float cmax[3] = { 1.000f, 0.933f, 0.900f }; /* Color for max positive spots. */
    /* Interpolate between them with {z}: */
    frgb_t clr;
    int c;
    if (fabs(z) <= 1.0e-5)
      { /* Map to center gray: */
        for (c = 0; c < 3; c++) { clr.c[c] = 0.500; }
      }
    else
      { /* Stretch argument away from 0: */
        double t = sqrt(fabs(z));
        /* Interpolate between {cmin} and {cmax} with {t} as ratio: */
        double s = 1 - t;
        for (c = 0; c < 3; c++) 
          { /* DeCasteljau's quadratic interpolation algorithm: */
            float v0t = (float)(s*cmin[c] + t*cmed[c]);
            float vt1 = (float)(s*cmed[c] + t*cmax[c]);
            float vtt = (float)(s*v0t + t*vt1);
            /* If {z} is negative, complement with respect to center  gray: */
            if (z < 0) { vtt = 1.000f - vtt; }
            clr.c[c] = vtt;
          }
      }
    return clr;
  }

frgb_t frgb_path_map_signed_1(double z, int cycles)
  { /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Compute {z} point along a spiral path: */
    frgb_t clr;
    int c;
    if (fabs(z) == 0.0)
      { /* Map to center gray: */
        for (c = 0; c < 3; c++) { clr.c[c] = 0.500; }
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

        for (c = 0; c < 3; c++) 
          { double Ac = qwht*1.0 + qell*(m[c] + au*u[c] + av*v[c]);
            double Bc = pm*m[c] + puv*u[c] + puv*v[c];
            double Cc = kspr*Ac + klin*Bc;
            /* If {z} is negative, complement with respect to center  gray: */
            if (z < 0) { Cc = 1.000 - Cc; }
            clr.c[c] = (float)Cc;
          }
      }
    
    return clr;
  }

frgb_t frgb_path_map_signed_2(double z, int cycles)
  { /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Compute {z} point along a spiral path: */
    frgb_t clr;
    int c;
    if (fabs(z) == 0.0)
      { /* Map to center gray: */
        for (c = 0; c < 3; c++) { clr.c[c] = 0.500; }
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

        for (c = 0; c < 3; c++) 
          { double Ac = qwht*1.0 + qell*(m[c] + au*u[c] + av*v[c]);
            double Cc = kspr*Ac + kmid*m[c];
            /* If {z} is negative, complement with respect to center  gray: */
            if (z < 0) { Cc = 1.000 - Cc; }
            clr.c[c] = (float)Cc;
          }
      }
    
    return clr;
  }
    
frgb_t frgb_path_map_signed(double z, int cycles, int style)
  { 
    switch(style) 
      {
        case 0: return frgb_path_map_signed_0(z, cycles); break;
        case 1: return frgb_path_map_signed_1(z, cycles); break;
        case 2: return frgb_path_map_signed_1(z, cycles); break;
        default: demand(FALSE, "invalid style"); 
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
