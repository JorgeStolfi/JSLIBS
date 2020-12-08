/* See pscolor.h */
/* Last edited on 2009-08-25 13:28:13 by stolfi */

#include <pswr_color.h>
#include <pswr_iso.h>

#include <affirm.h>

#include <math.h>
#include <stdlib.h>
#include <limits.h>

#define M_SQRT3 (1.73205080756887729352)
  /* Precomputed sqrt(3). */

void pswr_make_color_table
  ( double vStart,
    double vStep,
    int kMin,
    int kMax, 
    double RMin, double GMin, double BMin,
    double RZer, double GZer, double BZer,
    double RMax, double GMax, double BMax,
    int *NP, 
    double **RP, 
    double **GP,
    double **BP
  )
  {
    /* Create gradual color table: */
    int N = kMax - kMin + 2;     /* Number of color bands: */
    double *R = (double *)malloc(N*sizeof(double));
    double *G = (double *)malloc(N*sizeof(double));
    double *B = (double *)malloc(N*sizeof(double));
    double vMin = vStart + (kMin - 0.5)*vStep;
    double vMax = vStart + (kMax + 0.5)*vStep;
    
    int i;
    for (i = 0; i < N; i++)
      { int k = kMin + i;  /* Index of upper isoline of band {i}. */
        
        /* Compute relative function value {r} at mid-band, in [-1 _ +1]: */
        double r; 
        if (i == 0) 
          { r = -1.0; }
        else if (i == N-1) 
          { r = +1.0; }
        else
          { /* Function values at lower and upper isolines of band: */
            double vLo = pswr_level(vStart, vStep, k-1);
            double vHi = pswr_level(vStart, vStep, k);
            /* Compute nominal function value {vMd} at mid-band: */
            double vMd = (vLo + vHi)/2;
            affirm ((vMin < vMd) && (vMd < vMax), "bug");
            r = (vMd > 0.0 ? vMd/vMax : -vMd/vMin);
          }
          
        /* Intepolate between colors in equal perceptual distance: */
        double *Ri = &(R[i]), *Gi = &(G[i]), *Bi = &(B[i]);
        if (r > 0.0)
          { pswr_interpolate_colors(+r, RZer,GZer,BZer, RMax,GMax,BMax, Ri,Gi,Bi); }
        else
          { pswr_interpolate_colors(-r, RZer,GZer,BZer, RMin,GMin,BMin, Ri,Gi,Bi); }
      }
      
    *NP = N;
    *RP = R; *GP = G; *BP = B;
  }

void pswr_interpolate_colors
  ( double r, 
    double R0, double G0, double B0,
    double R1, double G1, double B1,
    double *R, double *G, double *B
  )
  {
    /* Compute brightnesses {Y0,Y1} of {R0,G0,B0} and {R1,G1,B1}: */ 
    double Y0 = 0.299*R0 + 0.587*G0 + 0.114*B0;
    double Y1 = 0.299*R1 + 0.587*G1 + 0.114*B1;
    /* Interpolate brightness {Y} in cubic-root scale: */
    double w = (1-r)*cbrt(Y0) + r*cbrt(Y1);
    double Y = w*w*w;
    if (Y < 0) { Y = 0; }
    if (Y > 1) { Y = 1; }
    
    /* Compute limiting chroma vectors: */ 
    double dR0 = R0 - Y0, dG0 = G0 - Y0, dB0 = B0 - Y0;
    double dR1 = R1 - Y1, dG1 = G1 - Y1, dB1 = B1 - Y1;
    /* Interpolate chroma vector {dR,dG,dB} in linear scale: */
    double dR = (1-r)*dR0 + r*dR1;
    double dG = (1-r)*dG0 + r*dG1;
    double dB = (1-r)*dB0 + r*dB1;
    
    /* Compute max {h <= 1} h.t. {(Y,Y,Y) + h*(dR,dG,dB)} is ok: */
    double h = 1.0, lim;
    lim = (dR < 0 ? 0 : 1) - Y; if (h*dR*dR > lim*dR) { h = lim/dR; }
    lim = (dG < 0 ? 0 : 1) - Y; if (h*dG*dG > lim*dG) { h = lim/dG; }
    lim = (dB < 0 ? 0 : 1) - Y; if (h*dB*dB > lim*dB) { h = lim/dB; }
    /* Now add the clipped chroma component: */
    (*R) = Y + h*dR;
    (*G) = Y + h*dG;
    (*B) = Y + h*dB;
  }

void pswr_color_scale_1
  ( double fs, 
    double Rs, double Gs, double Bs,
    double Y0,
    double *R, double *G, double *B
  )
  {
    /* Compute brightness {Ys} of {Rs,Gs,Bs}: */ 
    double Ys = 0.299*Rs + 0.587*Gs + 0.114*Bs;
    /* Compute the pseudo-saturation {m} of {fs} relative to the unit 1-ball: */ 
    double m = fabs(fs);
    /* Compute brightness {Y} of target color: */ 
    double bias = 0.125;
    double r = (Ys + bias)/(Y0 + bias);
    double Y = (Y0 + bias)*exp(m*log(r)) - bias;
    if (Y < 0) { Y = 0; }
    if (Y > 1) { Y = 1; }
    /* Compute limiting chroma vector (for {m = 1}): */ 
    double Rd = Rs - Ys, Gd = Gs - Ys, Bd = Bs - Ys;
    if (fs < 0) { Rd = -Rd; Gd = -Gd; Bd = -Bd; }
    /* Compute max {h <= 1} such that {(Y,Y,Y) + h*(Rd,Gd,Bd)} is ok: */
    double h = 1.0, lim;
    lim = (Rd < 0 ? 0 : 1) - Y; if (h*Rd*Rd > lim*Rd) { h = lim/Rd; }
    lim = (Gd < 0 ? 0 : 1) - Y; if (h*Gd*Gd > lim*Gd) { h = lim/Gd; }
    lim = (Bd < 0 ? 0 : 1) - Y; if (h*Bd*Bd > lim*Bd) { h = lim/Bd; }
    /* Now add the chroma component to {(Y,Y,Y)}: */
    m *= h;
    (*R) = Y + m*Rd;
    (*G) = Y + m*Gd;
    (*B) = Y + m*Bd;
  }

void pswr_color_scale_2
  ( double fs, double ft, 
    double Rs, double Gs, double Bs,
    double Y0,
    double *R, double *G, double *B
  )
  {
    /* Convert {Rs,Gs,Bs} to YUV coordinates: */ 
    double Ys = 0.299*Rs + 0.587*Gs + 0.114*Bs;
    double Us = 0.284*(Rs - Gs);
    double Vs = - 0.139*Rs - 0.085*Gs + 0.224*Bs;
    /* Compute the zero-brightness component {Rp,Gp,Bp} of {Rs,Gs,Bs}: */ 
    double Rp = Rs - Ys, Gp = Gs - Ys, Bp = Bs - Ys;
    /* Compute a zero-brightness vector {Rq,Gq,Bq} orthogonal to {Rp,Gp,Bp}: */ 
    double Ut = -Vs, Vt = Us;
    double Rq = + 2.2214430220*Ut - 0.50943800940*Vt;
    double Gq = - 1.2996837380*Ut - 0.50943800940*Vt;
    double Bq = + 0.8853011714*Ut + 3.95484770500*Vt;
    /* Compute the pseudo-saturation {m} of {fs} relative to the unit 2-ball: */ 
    double m = hypot(fs, ft);
    /* Compute the brightness {Y} of the target color: */ 
    double bias = 0.125;
    double r = (Ys + bias)/(Y0 + bias);
    double Y = (Y0 + bias)*exp(m*log(r)) - bias;
    if (Y < 0) { Y = 0; }
    if (Y > 1) { Y = 1; }
    /* Compute the limiting chroma vector (for {m = 1}): */
    double Rd = fs/m * Rp + ft/m * Rq;
    double Gd = fs/m * Gp + ft/m * Gq;
    double Bd = fs/m * Bp + ft/m * Bq;
    /* Compute the maximum {h <= 1} s.t. {(Y,Y,Y) + h*(Rd,Gd,Bd)} is ok: */
    double h = 1.0, lim;
    lim = (Rd < 0 ? 0 : 1) - Y; if (h*Rd*Rd > lim*Rd) { h = lim/Rd; }
    lim = (Gd < 0 ? 0 : 1) - Y; if (h*Gd*Gd > lim*Gd) { h = lim/Gd; }
    lim = (Bd < 0 ? 0 : 1) - Y; if (h*Bd*Bd > lim*Bd) { h = lim/Bd; }
    /* Now add the chroma component to {(Y,Y,Y)}: */
    m *= h;
    (*R) = Y + m*Rd;
    (*G) = Y + m*Gd;
    (*B) = Y + m*Bd;
  }
  
void pswr_color_scale_3
  ( double fs, double ft, double fu, 
    double Ymin,
    double Ymax,
    double *R, double *G, double *B
  )
  {
    (*R) = (1 - fs)/2*Ymin + (fs + 1)/2*Ymax;
    (*G) = (1 - ft)/2*Ymin + (ft + 1)/2*Ymax;
    (*B) = (1 - fu)/2*Ymin + (fu + 1)/2*Ymax;
  }

