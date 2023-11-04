/* See {frgb_path.h}. */
/* Last edited on 2023-03-06 19:42:45 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <frgb.h>
#include <jsmath.h>
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

void frgb_path_choose_unsigned_0_parms(int32_t cycles, int32_t *fR_P, double *aR_P, int32_t *fB_P, double *aB_P);
  /* Chooses the parameters {fR,aR,fB,aB} for the {frgb_path_map_unsigned_0} function,
    based on the {cycles} parameter. */

double frgb_path_clip_unsigned(double z);
  /* Clips {z} to the range {[0_1]}. */

double frgb_path_clip_signed(double z);
  /* Clips {z} to the range {[-1_+1]}. */

#define YR  (0.298911)
#define YG  (0.586611)
#define YB  (0.114478)
  /* Luminance weights according to European TV standard. */

frgb_t frgb_path_map_unsigned_0(double z, int32_t cycles)
  { 
    /* This path is based on an ad-hoc formula that produces a smooth
      (C infinity)sequence of colors, with fairly varied hues and
      monotonically increasing {Y}.
      
      For each {z}, {R} and {B} channels are computed by ramp plus sinusoidal functions
      {a*T(z)*(1-cos(f*cycles*PI*z))/2 + z} with different (half-)frequencies {f=fR,fB},
      and different parameters {a=aR,aB}, where {T(z)} is suitable hump function 
      that is zero at {z=0} and {z=1}.
     
      The {G} channels is then computed so that the brightness {Y} is
      equal to {z}: namely, {G = (z - YR*R - YB*B)/YG} where
      {YR,YG,YB} are the {Y} values of pure red, pure green, and pure
      blue.
       
      The parameters {aR} and {aB} are chosen so that the 
      resulting triple {R,G,B} is always inside the unit cube. Thus 
      there is no need for clipping, and the path remains C-infinity. */

    /* Reduce the path parameter to {[0_1]}: */
    z = frgb_path_clip_unsigned(z);

    /* Choose the half-period counts {fR,fB} and : */
    int32_t fR,fB;
    double aR, aB;
    frgb_path_choose_unsigned_0_parms(cycles, &fR,&aR,&fB,&aB);

    /* Compute the amplitude {T}: */
    double S = 0.5*(1 - cos(2*M_PI*z));
    double T = S*(2 - S);

    /* Compute {R} and {B}: */
    double R = z + aR*T*(1 - cos(fR*M_PI*z))/2;
    double B = z + aB*T*(1 - cos(fB*M_PI*z))/2;
    
    /* Compute {G} so that the brightness (European TV {Y}) is {z}: */
    double G = 1.0e-8 + (1.0 - 2.0e-8)*(z - YR*R - YB*B)/YG;
    fprintf(stderr, "%14.8f  %14.8f %14.8f %14.8f\n", z, R, G, B);
    
    assert((R >= 0.0) && (R <= 1.0));
    assert((G >= 0.0) && (G <= 1.0));
    assert((B >= 0.0) && (B <= 1.0));
    
    return (frgb_t){{ (float)R, (float)G, (float)B }};
  }

void frgb_path_choose_unsigned_0_parms(int32_t cycles, int32_t *fR_P, double *aR_P, int32_t *fB_P, double *aB_P)
  { int32_t fR= -1, fB = -1;
    double aR,aB;
    if (cycles == 0)
      { fR = 0; aR = 0; fB = 0; aB = 0; }
    else if (cycles < 0)
      { switch(-cycles) 
          { case 1: 
              fR = 1; aR = 0.383; fB = 3; aB = 0.663; break;
            case 2: 
              fR = 3; aR = 0.661; fB = 5; aB = 0.120; break;
            case 3: 
              fR = 5; aR = 0.396; fB = 7; aB = 0.336; break;
            case 4: 
              fR = 5; aR = 0.396; fB = 9; aB = 0.338; break;
            case 5: 
              fR = 5; aR = 0.396; fB = 9; aB = 0.338; break;
            default:
              break;
          }
      }
    else
      { switch(cycles) 
          { case 1: 
              fR =  3; aR = 0.655; fB = 1; aB = 0.335; break;
            case 2: 
              fR =  5; aR = 0.396; fB = 3; aB = 0.663; break;
            case 3: 
              fR =  7; aR = 0.336; fB = 5; aB = 0.396; break;
            case 4: 
              fR =  9; aR = 0.338; fB = 5; aB = 0.396; break;
            default:
              break;
          }
      }
      
    if (fR < 0)
      { /* Must do a general case: */
        fR = 2*abs(cycles) + 1;
        double phi = (sqrt(5)-1)/2;
        fB = (int32_t)floor(phi*fR);
        fB = fB + 1 - (fB%2);
        while ((fB == fR) || (fB <= 0) || (gcd(fR,fB) != 1)) { fB = fB + 2; }
        if (cycles < 0) { int32_t t = fR; fR = fB; fB = t; }
        /* Risk it: */
        aR = 0.330; aB = 0.330;
      }
      
    (*fR_P) = fR;
    (*aR_P) = aR;
    (*fB_P) = fB;
    (*aB_P) = aB;
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

frgb_t frgb_path_map_signed_0(double z)
  { /* A simple S-shaped path in the cube.  Ignores the {cycles} argument. */
    /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Choose extrmal colors: */
    float cmin[3] = { 0.500f, 0.500f, 0.500f }; /* Middle gray. */
    float cmed[3] = { 1.000f, 0.333f, 0.000f }; /* Bézier middle color. */
    float cmax[3] = { 1.000f, 0.933f, 0.900f }; /* Color for max positive spots. */
    /* Interpolate between them with {z}: */
    frgb_t clr;
    int32_t c;
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

frgb_t frgb_path_map_signed_1(double z, int32_t cycles)
  { /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Compute {z} point along a spiral path: */
    frgb_t clr;
    int32_t c;
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

        /* Map {abs(z)} to a point {au,av} on the {[-1_+1]×[0_+1]} square: */
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

frgb_t frgb_path_map_signed_2(double z, int32_t cycles)
  { /* Reduce the path parameter to {[-1_+1]}: */
    z = frgb_path_clip_signed(z);
    /* Compute {z} point along a spiral path: */
    frgb_t clr;
    int32_t c;
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

        /* Map {abs(z)} to a point {au,av} on the {[-1_+1]×[-1_+1]} square: */
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
    
void frgb_path_polygon_corners(double z, int32_t *nf_P, frgb_t f[])
  { 
    /* Check {z} range: */
    demand((z >= 0) & (z <= 1.0), "invalid brightness value");
    
    /* Intersect the plane {Y = z} with the edges of the RGB unit cube: */
    int32_t nf = 0;
    if (z < 0) 
      { /* No intersection. */ }
    else if (z == 0)
      { /* Just the black (K) corner: */
        f[nf] = (frgb_t){{ 0, 0, 0 }}; nf++;
      }
    else if (z <= YB) /* 0.0~0.1 */
      { /* Triangle; edges are K-R,K-G,K-B: */
        f[nf] = (frgb_t){{ (float)(z/YR), 0, 0 }}; nf++;
        f[nf] = (frgb_t){{ 0, (float)(z/YG), 0 }}; nf++;
        f[nf] = (frgb_t){{ 0, 0, (float)(z/YB) }}; nf++;
      }
    else if (z <= YR) /* 0.1~0.3 */
      { /* Quadrilateral; edges are K-R,K-G,B-C,B-M: */
        f[nf] = (frgb_t){{ (float)(z/YR), 0, 0 }};      nf++;
        f[nf] = (frgb_t){{ 0, (float)(z/YG), 0 }};      nf++;
        f[nf] = (frgb_t){{ 0, (float)((z-YB)/YG), 1 }}; nf++;
        f[nf] = (frgb_t){{ (float)((z-YB)/YR), 0, 1 }}; nf++;
      }
    else if (z < YR+YB) /* 0.3~0.4 */
      { /* Pentagon; edges are R-Y,K-G,B-C,B-M,R-M: */
        f[nf] = (frgb_t){{ 1, (float)((z-YR)/YG), 0 }}; nf++;
        f[nf] = (frgb_t){{ 0, (float)(z/YG), 0 }};      nf++;
        f[nf] = (frgb_t){{ 0, (float)((z-YB)/YG), 1 }}; nf++;
        f[nf] = (frgb_t){{ (float)((z-YB)/YR), 0, 1 }}; nf++;
        f[nf] = (frgb_t){{ 1, 0, (float)((z-YR)/YB) }}; nf++;
      }
    else if (z <= YG) /* 0.4~0.6 */
      { /* Quadrilateral; edges are R-Y,K-G,B-C,M-W: */
        f[nf] = (frgb_t){{ 1, (float)((z-YR)/YG), 0 }};    nf++;
        f[nf] = (frgb_t){{ 0, (float)(z/YG), 0 }};         nf++;
        f[nf] = (frgb_t){{ 0, (float)((z-YB)/YG), 1 }};    nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR-YB)/YG), 1 }}; nf++;
      }
    else if (z < YG+YB) /* 0.6~0.7 */
      { /* Pentago; edges are R-Y,G-Y,G-C,B-C,M-W: */
        f[nf] = (frgb_t){{ 1, (float)((z-YR)/YG), 0 }};    nf++;
        f[nf] = (frgb_t){{ (float)((z-YG)/YR), 1, 0 }};    nf++;
        f[nf] = (frgb_t){{ 0, 1, (float)((z-YG)/YB) }};    nf++;
        f[nf] = (frgb_t){{ 0, (float)((z-YB)/YG), 1 }};    nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR-YB)/YG), 1 }}; nf++;
      }
    else if (z < YR+YG) /* 0.7~0.9 */
      { /* Quadrilateral; edges are R-Y,G-Y,C-W,M-W: */
        f[nf] = (frgb_t){{ (float)((z-YG)/YR), 1, 0 }};    nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR)/YG), 0 }};    nf++;
        f[nf] = (frgb_t){{ (float)((z-YG-YB)/YR), 1, 1 }}; nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR-YB)/YG), 1 }}; nf++;
      }
    else if (z < 1.0) /* 0.9~1.0 */
      { /* Triangle; edges are Y-W,C-W,M-W: */
        f[nf] = (frgb_t){{ (float)((z-YG-YB)/YR), 1, 1 }}; nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR-YB)/YG), 1 }}; nf++;
        f[nf] = (frgb_t){{ 1, 1, (float)((z-YG-YR)/YB) }}; nf++;
      }
    else if (z == 1.0)
      { /* Single point, the white (W) corner: */
        f[nf] = (frgb_t){{ 1, 1, 1 }}; nf++;
      }
    else
      { /* No intersection. */ }
      
    if (nf >= 2)
      { /* Eliminate consecutive corners that are repeated by roundoff: */
        int32_t nf_new = 0; /* Number of distinct corners. */
        f[nf_new] = f[0]; nf_new++;
        for (int32_t k = 1; k < nf; k++)
          { if (! frgb_eq(&(f[k]), &(f[nf_new-1])))
              { f[nf_new] = f[k]; nf_new++; }
          }
        if (nf_new == 2)
          { /* Must be distinct only by roundoff: */
            nf_new = 1;
          }
        nf = nf_new;
      }
        
    (*nf_P) = nf;
  }

frgb_t frgb_path_map_signed(double z, int32_t cycles, int32_t style)
  { 
    switch(style) 
      {
        case 0: 
          demand(cycles == 1, "signed path style 0 does not support {cycles}"); 
          return frgb_path_map_signed_0(z);
          break;
        case 1: 
          return frgb_path_map_signed_1(z, cycles); 
          break;
        case 2: 
          return frgb_path_map_signed_1(z, cycles);
          break;
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

