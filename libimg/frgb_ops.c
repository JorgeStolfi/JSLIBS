/* See {frgb_ops.h}. */
/* Last edited on 2023-03-07 14:14:53 by stolfi */

/* Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#define _GNU_SOURCE
#include <math.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <frgb_ops.h>
#include <sample_conv.h>

#include <argparser.h>
#include <frgb.h>
#include <bool.h>
#include <fget.h>

int32_t frgb_DEBUG = 0;

frgb_t frgb_scale(double s, frgb_t *a)
  { frgb_t p;
    for (int32_t i = 0; i < 3; i++) { p.c[i] = (float)(s * a->c[i]); }
    return p;
  }

frgb_t frgb_shift(double d, frgb_t *a)
  { frgb_t p;
    for (int32_t i = 0; i < 3; i++) { p.c[i] = (float)(a->c[i] + d); }
    return p;
  }

frgb_t frgb_add(frgb_t *a, frgb_t *b)
  { frgb_t p;
    for (int32_t i = 0; i < 3; i++) { p.c[i] = (float)(a->c[i] + b->c[i]); }
    return p;
  }
  
frgb_t frgb_sub(frgb_t *a, frgb_t *b)
  { frgb_t p;
    for (int32_t i = 0; i < 3; i++) { p.c[i] = a->c[i] - b->c[i]; }
    return p;
  }
  
frgb_t frgb_mul(frgb_t *a, frgb_t *b)
  { frgb_t p;
    for (int32_t i = 0; i < 3; i++) { p.c[i] = a->c[i]*b->c[i]; }
    return p;
  }

frgb_t frgb_mix(double ca, frgb_t *a, double cb, frgb_t *b)
  { frgb_t p;
    for (int32_t i = 0; i < 3; i++) { p.c[i] = (float)(ca*a->c[i] + cb*b->c[i]); }
    return p;
  }

bool_t frgb_is_all_zeros(frgb_t *a)
  { return ((a->c[0] == 0) && (a->c[1] == 0) && (a->c[2] == 0)); }

bool_t frgb_is_all_ones(frgb_t *a)
  { return ((a->c[0] == 1) && (a->c[1] == 1) && (a->c[2] == 1)); }

bool_t frgb_eq(frgb_t *a, frgb_t *b)
  { return ((a->c[0] == b->c[0]) && (a->c[1] == b->c[1]) && (a->c[2] == b->c[2])); }

frgb_t frgb_parse(argparser_t *pp, double lo, double hi)
  { frgb_t p;
    for (int32_t i = 0; i < 3; i++)
      { p.c[i] = (float)argparser_get_next_double(pp, lo, hi); }
    return p;
  }

frgb_t frgb_parse_color(argparser_t *pp)
  { frgb_t p;
    double v[3];
    double scale;
    for (int32_t i = 0; i < 3; i++)
      { v[i] = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); }
    if (argparser_keyword_present_next(pp, "/"))
      { scale = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); }
    else
      { scale = 1.0; }
    for (int32_t i = 0; i < 3; i++)
      { p.c[i] = (float)(v[i]/scale); }
    return p;
  }

frgb_t frgb_read(FILE *rd, double lo, double hi)
  { frgb_t p;
    for (int32_t i = 0; i < 3; i++)
      { double pi = fget_double(rd);
        if ((pi < lo) || (pi > hi))
          { fprintf(stderr, "  component [%d] = %25.16e\n", i, pi);
            demand(FALSE, "bad {frgb_t} component value");
          }
        p.c[i] = (float)pi;
      }
    return p;
  }

frgb_t frgb_read_color(FILE *rd)
  { frgb_t p;
    double v[3];
    double scale;
    for (int32_t i = 0; i < 3; i++)
      { v[i] = fget_double(rd); }
    fget_skip_spaces(rd);
    if (fget_test_char(rd, '/'))
      { scale = fget_double(rd); }
    else
      { scale = 1.0; }
    for (int32_t i = 0; i < 3; i++)
      { p.c[i] = (float)(v[i]/scale); }
    return p;
  }

double frgb_floatize(int32_t ival, int32_t maxval, double zero, double scale)
  { 
    return (float)(((double)ival - zero)/scale);
  }
  
double frgb_log_scale_gray(double x)
  { if (x < VAL_EPS) { x = VAL_EPS; }
    return log(x);
  }

void frgb_log_scale(frgb_t *p, int32_t chns)
  { for (int32_t i = 0; i < chns; i++)
      { double xi = p->c[i];
        if (xi < VAL_EPS) { xi = VAL_EPS; }
        p->c[i] = (float)log(xi);
      }
  }

double frgb_clip_gray(double x)
  {  if (x < 0.0)
      { return 0.0; }
    else if (x > 1.0)
      { return 1.0; }
    else
      { return x; }
  }
  
double frgb_gamma_encoding_gray(double y, double gamma, double bias)
{ return (double)sample_conv_gamma((float)y, 1/gamma, bias); }

double frgb_gamma_decoding_gray(double y, double gamma, double bias)
  { return (double)sample_conv_gamma((float)y, gamma, bias); }

#define BT_BIAS (sample_conv_BT709_BIAS)
  /* !!! This should be a parameter of {frgb_correct_arg}. !!! */

frgb_t frgb_correct_arg(frgb_t *p, frgb_t *inGamma, int32_t gray)
  { frgb_t res = *p;
    if (inGamma != NULL)
      { for (int32_t i = 0; i < 3; i++)
          { res.c[i] = sample_conv_gamma(res.c[i], inGamma->c[i], BT_BIAS); }
      }
    if (gray)
      { res.c[0] = res.c[1] = res.c[2] = (float)frgb_luminance_CIE_XYZrec601_1(&res); }
    return res;
  }
  
void frgb_clip_rgb(frgb_t *p)
  { float lum = (float)frgb_luminance_CIE_XYZrec601_1(p); /* Luminance */
    float obs = 1.0f - lum;  /* "Obscurance" = complement of luminance */
    /* Check for extreme values: */
    if (lum <= 0.0)
      { *p = frgb_Black; }
    else if (obs <= 0.0)
      { *p = frgb_White; }
    else
      { frgb_t q = (frgb_t){{ lum, lum, lum }};
        return frgb_clip_rgb_towards(p, &q);
      }
  }
  
void frgb_clip_rgb_towards_grey(frgb_t *p)
  { frgb_t q = (frgb_t){{ 0.5f, 0.5f, 0.5f }};
    return frgb_clip_rgb_towards(p, &q);
  }

void frgb_clip_rgb_towards(frgb_t *p, frgb_t *q)
  { /* Compute the min {s >= 1} such that {q + (p - q)/s} is in the cube: */
    double s = 1.0;
    for (int32_t i = 0; i < 3; i++)
      { double pi = p->c[i];
        double qi = q->c[i];
        double qf = 1 - qi;
        demand((qi > 0) && (qf > 0), "q is not in the cube's interior");
        double di = pi - qi;
        if (s*qi < -di) 
          { s = -di/qi; }
        else if (s*qf < +di)
          { s = +di/qf; }
      }
    if (s > 1.0)
      { /* Pull {p} towards {q} until inside the unit cube: */
        for (int32_t i = 0; i < 3; i++) 
          { double ri = q->c[i] + (p->c[i] - q->c[i])/s;
            /* Guarding against roundoff: */
            if (ri < 0) { ri = 0; }
            if (ri > 1) { ri = 1; }
            p->c[i] = (float)ri;
          }
      }
  }

double frgb_apply_kappa_gray(double y, double kappa)
  { if (y <= 0.0)
      { return 0.0; }
    else if (y >= 1.0)
      { return 1.0; }
    else if (kappa == 1.0)
      { return y; }
    else
      { return kappa*y/((kappa-1.0)*y + 1.0); }
  }

void frgb_apply_glob_kappa_sat_clip(frgb_t *p, double kap, double satf)
  { double lumI = frgb_luminance_CIE_XYZrec601_1(p);
    double obsI = 1.0 - lumI;  /* "Obscurance" = complement of luminance */
    double satI = 0.0;
    double lumO, obsO;
    double satO = 0.0;
    /* Compute output luminance {lumO} and its complement {obsO}: */
    if (kap != 1.0)
      { lumO = frgb_apply_kappa_gray(lumI, kap); obsO = 1.0 - lumO; }
    else
      { lumO = lumI; obsO = obsI; }
    /* Check for extreme values: */
    if (lumO <= 0.0)
      { *p = frgb_Black; }
    else if (obsO <= 0.0)
      { *p = frgb_White; }
    else
      { /* Replace {p} by its chrominance component {lumI*White - p}
          If {satf < 1.0}, scale the chrominance right away,
          and set {satf} to 1.0; else leave this adjustment for later.
          Then compute the input saturation {satI}, relative to the maximum
          possible saturation for {p}'s hue and luminosity {lumI}:
            { satI = MAX{ s : lumI*White + p/s IN [0_1]^3 } }
          Then clip {satI} to the range {[0_1]}.
          Compute also the saturation {satO} of the same chrominance
          vector, relative to the maximum saturation for {p}'s hue and  
          for the desired luminosity {lumO}:
            { satO = MAX{ s : lumO*White + p/s IN [0_1]^3 } }
        */
        for (int32_t i = 0; i < 3; i++)
          { double pi = p->c[i];
            pi = pi - lumI;
            if (satf < 1.0) { pi = pi * satf; }
            if (satI < 1.0)
              { if ((pi > +obsI) || (pi < -lumI))
                  { satI = 1.0; }
                else if (satI*lumI < -pi)
                  { satI = -pi/lumI; }
                else if (satI*obsI < +pi)
                  { satI = +pi/obsI; }
              }
            if (satO*lumO < -pi)
              { satO = -pi/lumO; }
            else if (satO*obsO < +pi)
              { satO = +pi/obsO; }
            p->c[i] = (float)pi;
          }
        if (satf < 1.0) { satf = 1.0; }
        /* Compute output pixel {p} from new luminance {lumO} plus the 
          input chroma vector {p}, the latter scaled by
          {f == r/((1-satI) + r*satO)} where {r == satf*lumO/lumI}: */
        { double r = satf*lumO/lumI;
          if (r == 0.0)
            { *p = (frgb_t){{ (float)lumO, (float)lumO, (float)lumO }}; }
          else
            { double f = r/((1.0 - satI) + r*satO);
              for (int32_t i = 0; i < 3; i++) { p->c[i] = (float)(lumO + f * (double)(p->c[i])); }
            }
        }
      }
  }

int32_t frgb_quantize(double fval, double zero, double scale, int32_t maxval)
  { if (fval == +INF)
      { return maxval; }
    else if (fval == -INF)
      { return 0; }
    else
      { double ival = fval*scale + zero;
        if (isnan(ival))
          { return maxval/2; }
        else if (ival <= 0.0)
          { return 0; }
        else if (ival >= (double)maxval)
          { return maxval; }
        else
          { return (int32_t)(floor(ival + 0.5)); }
      }
  }
  
int32_t frgb_dequal(double *a, double *b, int32_t chns)
  { for (int32_t i = 0; i < chns; i++)
      { if (a[i] != b[i]) return 0; }
    return 1;
  }
  
int32_t frgb_fequal(float *a, float *b, int32_t chns)
  { for (int32_t i = 0; i < chns; i++)
      { if (a[i] != b[i]) return 0; }
    return 1;
  }

void frgb_print(FILE *f, char *pref, frgb_t *p, int32_t chns, char *fmt, char *suff)
  { demand(chns <= 3, "invalid channel count"); 
    fprintf(f, "%s", pref);
    for (int32_t k = 0; k < chns; k++) 
      { if (k > 0) { fputc(' ', f); }
        fprintf(f, fmt, p->c[k]);
      }
    fprintf(f, "%s", suff);
  }

void frgb_print_int_pixel(FILE *f, char *pref, int32_t *p, int32_t chns, char *suff)
  { fprintf(f, "%s", pref);
    for (int32_t k = 0; k < chns; k++) 
      { fprintf(f, "%s%5d", (k == 0 ? "" : " "), p[k]); }
    fprintf(f, "%s", suff);
  }

void frgb_debug(char *label, int32_t col, int32_t row, frgb_t *p, int32_t chns, char *tail)
  { 
    if (frgb_DEBUG) 
      { fprintf(stderr, "%s[%3d][%3d] = ", label, row, col);
        frgb_print(stderr, "( ", p, chns, "%10.7f", ")");
        fprintf(stderr, "%s", tail);
      }
  }

void frgb_debug_int_pixel(char *label, int32_t col, int32_t row, int32_t *p, int32_t chns, char *tail)
  { 
    if (frgb_DEBUG) 
      { fprintf(stderr, "%s[%3d][%3d] = ", label, row, col);
        frgb_print_int_pixel(stderr, "( ", p, chns, " )");
        fprintf(stderr, "%s", tail);
      }
  }

/* COLORSPACE CONVERSIONS */

/*
    Some RGB color space parameters (from Susstrunk, Buckley and Swen 2005)
    as reproduced in Wikipedia:

      Color Space        | WhtPt | Primaries
                         |       |  xR    | yR      | xG      | yG     | xB     | yB
      ISO RGB            | ANY   | 
      Extended ISO RGB   | ANY   | floating
      sRGB, HDTV[1]      | D65   | 0.64   | 0.33    | 0.30    | 0.60   | 0.15   | 0.06
      ROMM RGB           | D50   | 0.7347 | 0.2653  | 0.1596  | 0.8404 | 0.0366 | 0.0001
      Adobe RGB 98       | D65   | 0.64   | 0.34    | 0.21    | 0.71   | 0.15   | 0.06
      Apple RGB          | D65   | 0.625  | 0.34    | 0.28    | 0.595  | 0.155  | 0.07
      NTSC [2]           | C     | 0.67   | 0.33    | 0.21    | 0.71   | 0.14   | 0.08
      NTSC (1979) [3]    | D65   | 0.63   | 0.34    | 0.31    | 0.595  | 0.155  | 0.07
      PAL/SECAM [4]      | D65   | 0.64   | 0.33    | 0.29    | 0.60   | 0.15   | 0.06

      [1] ITU-R BT.709-3
      [2] FCC 1953 
      [3] SMPTE C, SMPTE-RP 145
      [4] EBU 3213, ITU-R BT.470-6
      
    It is not clear which RGB space is assumed by the formulas below.
    Someday we should replace this interface by one with the CIE XYZ
    system as the central standard, and the various RGB systems as 
    peripheral.
*/

void frgb_to_CIE_XYZrec601_1(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    
    double X = frgb_XR*R + frgb_XG*G + frgb_XB*B;
    double Y = frgb_YR*R + frgb_YG*G + frgb_YB*B;
    double Z = frgb_ZR*R + frgb_ZG*G + frgb_ZB*B;
    
    p->c[0] = (float)X;
    p->c[1] = (float)Y; 
    p->c[2] = (float)Z;
  }
  
void frgb_from_CIE_XYZrec601_1(frgb_t *p)
  { 
    double X = (double)p->c[0];
    double Y = (double)p->c[1];
    double Z = (double)p->c[2];

    double R = +1.910027*X -0.532464*Y -0.288214*Z;
    double G = -0.984645*X +1.999130*Y -0.028308*Z;
    double B = +0.058309*X -0.118385*Y +0.897608*Z;

    p->c[0] = (float)R;
    p->c[1] = (float)G; 
    p->c[2] = (float)B;
  }

double frgb_luminance_CIE_XYZrec601_1(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    return frgb_YR*R + frgb_YG*G + frgb_YB*B;
  }

void frgb_to_CIE_XYZccir709(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    
    double X = +0.412411*R +0.357585*G +0.180454*B;
    double Y = +0.212649*R +0.715169*G +0.072182*B;
    double Z = +0.019332*R +0.119195*G +0.950390*B;
   
    p->c[0] = (float)X;
    p->c[1] = (float)Y; 
    p->c[2] = (float)Z;
  }
  
void frgb_from_CIE_XYZccir709(frgb_t *p)
  { 
    double X = (double)p->c[0];
    double Y = (double)p->c[1];
    double Z = (double)p->c[2];

    double R = +3.240811*X -1.537310*Y -0.498586*Z;
    double G = -0.969241*X +1.875966*Y +0.041554*Z;
    double B = +0.055638*X -0.204007*Y +1.057130*Z;

    p->c[0] = (float)R;
    p->c[1] = (float)G; 
    p->c[2] = (float)B;
  }

double frgb_luminance_CIE_XYZccir709(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    return +0.212649*R +0.715169*G +0.072182*B;
  }


void frgb_to_CIE_XYZitu_D65(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    
    double X = +0.430574*R +0.341550*G +0.178325*B;
    double Y = +0.222015*R +0.706655*G +0.071330*B;
    double Z = +0.020183*R +0.129553*G +0.939180*B;  
   
    p->c[0] = (float)X;
    p->c[1] = (float)Y; 
    p->c[2] = (float)Z;
  }
  
void frgb_from_CIE_XYZitu_D65(frgb_t *p)
  { 
    double X = (double)p->c[0];
    double Y = (double)p->c[1];
    double Z = (double)p->c[2];

    double R = +3.063219*X -1.393326*Y -0.475801*Z;
    double G = -0.969245*X +1.875968*Y +0.041555*Z;
    double B = +0.067872*X -0.228833*Y +1.069251*Z;

    p->c[0] = (float)R;
    p->c[1] = (float)G; 
    p->c[2] = (float)B;
  }

double frgb_luminance_CIE_XYZitu_D65(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    return +0.222015*R +0.706655*G +0.071330*B;
  }

void frgb_to_YUV(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    
    double Y = frgb_YR*R + frgb_YG*G + frgb_YB*B;
    double U = -0.147000*R - 0.289000*G +0.436000*B;
    double V = +0.615000*R - 0.515000*G -0.100000*B;
   
    p->c[0] = (float)Y;
    p->c[1] = (float)U; 
    p->c[2] = (float)V;
  }
  
void frgb_from_YUV(frgb_t *p)
  { 
    double Y = (double)p->c[0];
    double U = (double)p->c[1];
    double V = (double)p->c[2];

    double R = +1.000000*Y -0.001164*U +1.139704*V;
    double G = +1.000000*Y -0.395735*U -0.580624*V;
    double B = +1.000000*Y +2.030875*U -0.000606*V;

    p->c[0] = (float)R;
    p->c[1] = (float)G; 
    p->c[2] = (float)B;
  }

double frgb_get_Y(frgb_t *p)
  { double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    
    double Y = frgb_YR*R + frgb_YG*G + frgb_YB*B;
    
    return Y;
  }

double frgb_get_Y_pbm(frgb_t *p)
  { double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    
    double Y = +0.299*R +0.587*G +0.114*B;
    
    return Y;
  }

void frgb_to_YIQ(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];

    double Y = frgb_YR*R + frgb_YG*G + frgb_YB*B;
    double I = +0.596018*R -0.274147*G -0.321871*B;
    double Q = +0.211646*R -0.522976*G +0.311330*B;
   
    p->c[0] = (float)Y;
    p->c[1] = (float)I; 
    p->c[2] = (float)Q;
  }
  
void frgb_from_YIQ(frgb_t *p)
  { 
    double Y = (double)p->c[0];
    double I = (double)p->c[1];
    double Q = (double)p->c[2];

    double R = +1.000000*Y +0.955920*I +0.620580*Q;
    double G = +1.000000*Y -0.271330*I -0.648223*Q;
    double B = +1.000000*Y -1.105629*I +1.701256*Q;
   
    p->c[0] = (float)R;
    p->c[1] = (float)G; 
    p->c[2] = (float)B;
  }


void frgb_to_YCbCr_601_1(frgb_t *p)
  { 
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];

    double Y  = frgb_YR*R + frgb_YG*G + frgb_YB*B;
    double Cb = -0.168777*R -0.331223*G +0.500000*B;
    double Cr = +0.500000*R -0.418357*G -0.081643*B;

    p->c[0] = (float)Y;
    p->c[1] = (float)Cb; 
    p->c[2] = (float)Cr;
  }
  
void frgb_from_YCbCr_601_1(frgb_t *p)
  { 
    double Y  = (double)p->c[0];
    double Cb = (double)p->c[1];
    double Cr = (double)p->c[2];

    double R = +1.000000*Y              +1.402178*Cr;
    double G = +1.000000*Y -0.345622*Cb -0.714488*Cr;
    double B = +1.000000*Y +1.771044*Cb;
   
    p->c[0] = (float)R;
    p->c[1] = (float)G; 
    p->c[2] = (float)B;
  }


void frgb_to_YUV_a(frgb_t *p)
  { 
    /* Old frgb_yuv conversion - matrix:
        [ 306,  601, 117]
        [ 291, -291,   0] / 1024
        [-142,  -87, 229]
      The result ranges in
        [0..1], [-291..+291]/1024, [-229..+229]/1024.
    */
    
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    
    double Y = frgb_YR*R + frgb_YG*G + frgb_YB*B;
    double U = +0.132812*R -0.132812*G;
    double V = -0.165039*R +0.038086*G +0.126953*B;

    p->c[0] = (float)Y;
    p->c[1] = (float)U; 
    p->c[2] = (float)V;
  }
  
void frgb_from_YUV_a(frgb_t *p)
  { 
    double Y = (double)p->c[0];
    double U = (double)p->c[1];
    double V = (double)p->c[2];

    double R = +1.000000*Y +4.158265*U -0.901735*V;
    double G = +1.000000*Y -3.371175*U -0.901735*V;
    double B = +1.000000*Y +6.417103*U +6.975196*V;

    p->c[0] = (float)R;
    p->c[1] = (float)G; 
    p->c[2] = (float)B;
  }

void frgb_to_YUV_b(frgb_t *p)
  { 
    /* Conversion used in old ppmytouvx */
    
    double R = (double)p->c[0];
    double G = (double)p->c[1];
    double B = (double)p->c[2];
    
    double Y = frgb_YR*R + frgb_YG*G + frgb_YB*B;
    double U = +0.284000*R -0.284000*G;
    double V = -0.139000*R -0.085000*G +0.224000*B;

    p->c[0] = (float)Y;
    p->c[1] = (float)U; 
    p->c[2] = (float)V;
  }
  
void frgb_from_YUV_b(frgb_t *p)
  { 
    double Y = (double)p->c[0];
    double U = (double)p->c[1];
    double V = (double)p->c[2];

    double R = +1.000000*Y +2.218491*U -0.511063*V;
    double G = +1.000000*Y -1.302636*U -0.511063*V;
    double B = +1.000000*Y +0.882349*U +3.953223*V;

    p->c[0] = (float)R;
    p->c[1] = (float)G; 
    p->c[2] = (float)B;
  }

void frgb_YUV_to_yuv(frgb_t *p, double ybias)
  { 
    /* Scale {Y,U,V} by {s = (1 + ybias)/(Y + ybias)}. */
    double Y = p->c[0];
    double num = 1 + ybias;
    double den = Y + ybias;
    /* Avoid excessive scale: */
    double tiny = ybias/2;
    if (den < tiny) den = tiny;
    double scale = num/den;
    /* Apply scale to {Y,U,V}: */
    p->c[0] = (float)((double)(p->c[0])*scale);
    p->c[1] = (float)((double)(p->c[1])*scale);
    p->c[2] = (float)((double)(p->c[2])*scale);
  }

void frgb_YUV_from_yuv(frgb_t *p, double ybias)
  { 
    /* Unscale {y,u,v} by {s = (1 + ybias)/(Y + ybias)} where {Y = y/s}. */
    double y = p->c[0];
    double num = 1 - y + ybias;
    double den = ybias;
    /* avoid scale smaller than 1: */
    if (num < den) num = den;
    double scale = num/den;
    /* Apply scale to {Y,U,V}: */
    p->c[0] = (float)((double)(p->c[0])/scale);
    p->c[1] = (float)((double)(p->c[1])/scale);
    p->c[2] = (float)((double)(p->c[2])/scale);
  }

void frgb_YUV_to_Yuv(frgb_t *p, double ybias)
  { 
    /* Scale {Y,U,V} by {s = (1 + ybias)/(Y + ybias)}. */
    double Y = p->c[0];
    double num = 1 + ybias;
    double den = Y + ybias;
    /* Avoid excessive scale: */
    double tiny = ybias/2;
    if (den < tiny) den = tiny;
    double scale = num/den;
    /* Apply scale to {U,V} only: */
    p->c[1] = (float)((double)(p->c[1])*scale);
    p->c[2] = (float)((double)(p->c[2])*scale);
  }

void frgb_YUV_from_Yuv(frgb_t *p, double ybias)
  { 
    /* Unscale {y,u,v} by {s = (1 + ybias)/(Y + ybias)} where {Y = y/s}. */
    double y = p->c[0];
    double num = 1 - y + ybias;
    double den = ybias;
    /* avoid scale smaller than 1: */
    if (num < den) num = den;
    double scale = num/den;
    /* Apply scale to {U,V} only: */
    p->c[1] = (float)((double)(p->c[1])/scale);
    p->c[2] = (float)((double)(p->c[2])/scale);
  }

void frgb_to_HSV_CG(frgb_t *p)
  { /* Grab the coordinates {R,G,B} of {p}: */
    double R = p->c[0], G = p->c[1], B = p->c[2];
    /* Find the max and min coordinates {cmin,cmax}: */
    double cmin = fmin(R, fmin(G, B));
    double cmax = fmax(R, fmax(G, B));
    /* Saturation and value are easy: */
    double S = cmax - cmin;
    double V = cmax;
    double H;
    if (S == 0)
      { H = 0; }
    else if (cmax == R)
      { if (cmin == B)
          { H = 0 + (G - cmin) / S; }
        else
          { H = 5 + (cmax - B) / S; }
      }
    else if (cmax == G)
      { if (cmin == R)
          { H = 2 + (B - cmin) / S; }
        else
          { H = 1 + (cmax - R) / S; }
      }
    else /* (cmax == B) */
      { if (cmin == G)
          { H = 4 + (R - cmin) / S; }
        else
          { H = 3 + (cmax - G) / S; }
      }
    H /= 6;
    (*p) = (frgb_t){{ (float)H, (float)S, (float)V }};
  }
  
void frgb_from_HSV_CG(frgb_t *p)
  { /* Grab the coordinates {H,S,V} of {p}: */
    double H = p->c[0], S = p->c[1], V = p->c[2];
    /* Reduce {H} to the range [0_1) modulo 1: */
    H = H - floor(H);
    assert((H >= 0) && (H < 1));
    /* Scale {H} to the range [0_6], split into integer {q} and fraction {f}: */
    H = H * 6;
    int32_t q = (int32_t)floor(H);
    double f = H - q;
    /* Snap to integer to account for roundoff in {1/6}, {2/6}, etc: */
    float sixth = ((float)1)/((float)6);
    double eps = 4*fabs(((double)sixth) - (1.0/6.0));
    if (fabs(f - 0) < eps) { f = 0; }
    if (fabs(f - 1) < eps) { f = 0; q = (q + 1) % 6; }
    /* Compute {R,G,B} of pure color with hue {H}: */
    double R, G, B;
    switch (q)
      { case 0:  R = 1;   G = f;   B = 0;   break;
        case 1:  R = 1-f; G = 1;   B = 0;   break;
        case 2:  R = 0;   G = 1;   B = f;   break;
        case 3:  R = 0;   G = 1-f; B = 1;   break;
        case 4:  R = f;   G = 0;   B = 1;   break;
        case 5:  R = 1;   G = 0;   B = 1-f; break;
        default: 
        /* Can't happen: */ 
        assert(FALSE);
      }
    /* Adjust by saturation and value: */
    R = V - (1-R)*S; 
    G = V - (1-G)*S; 
    B = V - (1-B)*S; 
    (*p) = (frgb_t){{ (float)R, (float)G, (float)B }};
  }

#define arg_UV_red (1.8054186283164418)
  /* Argument angle (radians) of the RGB red axis in the YUV system 
    projects on the UV plane, namely {U = -0.147, V = 0.615}. */

double frgb_H_from_UV(double U, double V)
  { double H;
    if ((U == 0) && (V == 0))
      { H = 0.0; }
    else
      { H = (atan2(V,U) - arg_UV_red)/(2*M_PI);
        while (H > 1) { H = H - 1; }
        while (H < 0) { H = H + 1; }
      }
    return H;
  }
  
void frgb_H_to_uv(double H, double *u, double *v)
  { double Hrad = H*2*M_PI + arg_UV_red;
    *u = cos(Hrad); *v = sin(Hrad);
  }

double frgb_get_H(frgb_t *p)
  { /* Convert to YUV and grab the U,V coordinates: */
    frgb_t q = (*p);
    frgb_to_YUV(&q);
    double U = q.c[1], V = q.c[2];
    /* The hue {H} is the argument of the {U,V} vector, scaled to period 1: */
    return frgb_H_from_UV(U, V);
  }

void frgb_YUV_to_YHS(frgb_t *p)
  { double Y = p->c[0];
    if ((Y <= 0) || (Y >= 1))
      { (*p) = (frgb_t){{ (float)Y, 0, 0 }}; }
    else
      { double U = p->c[1], V = p->c[2];
        double S = hypot(U, V)/Y;
        double H = frgb_H_from_UV(U, V);
        /* Fix possible roundoff overflows: */
        if (H >= 1.0) { H = H - 1; }
        if (H < 0.0) { H = H + 1; }
        (*p) = (frgb_t){{ (float)Y, (float)H, (float)S }};
      }
  }
 
void frgb_YHS_to_YUV(frgb_t *p)
  { 
    double Y = p->c[0];
    if ((Y <= 0) || (Y >= 1))
      { (*p) = (frgb_t){{ (float)Y, 0, 0 }}; }
    else
      { double S = p->c[2];
        if (S == 0)
          { (*p) = (frgb_t){{ (float)Y, 0, 0 }}; }
        else
          { demand(S > 0, "invalid saturation value"); 
            double H = p->c[1];
            double u, v;
            frgb_H_to_uv(H, &u, &v);
            (*p) = (frgb_t){{ (float)Y, (float)(Y*S*u), (float)(Y*S*v) }};
          }
      }
  }
  
double frgb_T_from_YUV(double Y, double U, double V)
  { if ((Y <= 0) || (Y >= 1))
      { return 0.0; }
    else
      { /* Compute a  zero-luminance RGB vector in direction {U,V}: */
        frgb_t Q = (frgb_t){{ 0.0, (float)U, (float)V  }};
        frgb_from_YUV(&Q);
        double R = Q.c[0], G = Q.c[1], B = Q.c[2];
        
        double T = 0;
        if (R != 0) { double TR = (R > 0 ? R/(1-Y) : -R/Y); if (TR > T) { T = TR; } }
        if (G != 0) { double TG = (G > 0 ? G/(1-Y) : -G/Y); if (TG > T) { T = TG; } }
        if (B != 0) { double TB = (B > 0 ? B/(1-Y) : -B/Y); if (TB > T) { T = TB; } }
        assert(isfinite(T) && (T >= 0));
        return T;
      }        
  }

void frgb_to_HTY(frgb_t *p)
  { /* Convert to YUV and grab those coordinates: */
    frgb_to_YUV(p);
    double Y = p->c[0], U = p->c[1], V = p->c[2];
    if ((Y <= 0) || (Y >= 1))
      { (*p) = (frgb_t){{ 0.0, 0.0, (float)Y }}; }
    else 
      { /* The hue {H} is the argument of the {U,V} vector, scaled to period 1: */
        double H = frgb_H_from_UV(U, V);
        assert((H >= 0) && (H <= 1.0));
        /* Compute the relative saturation {T}: */
        double T = frgb_T_from_YUV(Y, U, V);
        /* Repack: */
        (*p) = (frgb_t){{ (float)H, (float)T, (float)Y }};
      }
  }

void frgb_from_HTY(frgb_t *p)
  { /* Grab the coordinates {H,T,Y} of {p}: */
    double H = p->c[0], T = p->c[1], Y = p->c[2];
    demand(T >= 0, "invalid relative saturation {T}");
    if ((Y <= 0) || (Y >= 1.0) || (T == 0.0))
      { /* Gray color: */
        (*p) = (frgb_t) {{ (float)Y, (float)Y, (float)Y }};
      }
    else
      { /* Compute the {U,V} diretion {u,v} of hue {H}: */
        double u, v;
        frgb_H_to_uv(H, &u, &v);
        /* Compute the relative saturation of the unit {(u,v)} vector: */
        double t = frgb_T_from_YUV(Y, u, v);
        /* Compute the {(U,V)} coordinates from {t,T,u,v}: */
        assert(t > 0);
        double U = T*u/t, V = T*v/t;
        /* Repack as YUV: */
        (*p) = (frgb_t){{ (float)Y, (float)U, (float)V }}; 
        /* Convert to RGB: */
        frgb_from_YUV(p);
      }
  }
  
