/* Plots of bivariate functions with isolines or color bands. */
/* Last edited on 2009-08-24 22:31:54 by stolfi */

#ifndef pswr_iso_H
#define pswr_iso_H

#include <pswr.h>

/* 
  This module provides tools for painting a function defined in a
  region, with color banding or isolines (level curves). */
  
/* 
  LEVELS

  The isoline levels are an arithmetic progression specified by two
  parameters {vStart, vStep}: {v[k] = vStart + k*vStep}, with integer
  {k}. The step {vStep} must be positive. The integer {k} is called
  the /index/ of the isoline.
  
  Beware that, because of rondoff errors, the isoline spacing {v[k] -
  v[k-1]} may not be exactly uniform and equal to {vStep}. Beware also
  that the expression {vStart + k*vStep} may have different values
  depending on whether the compiler chooses to evaluate it in {double}
  or {double extended}.  */
  
double pswr_level(double vStart, double vStep, int k);
  /* The official level {v[k]} of the isoline with index {k}. */

int pswr_inf_isoline(double vStart, double vStep, double z);
  /* Returns the index of the first isoline at or below value {z},
    i.e. the largest integer {k} such that {v[k] \leq z}. */

int pswr_sup_isoline(double vStart, double vStep, double z);
  /* Returns the index of the first isoline at or above value {z},
    i.e. the smallest integer {k} such that {z \leq v[k]}. */

/* 
  BANDS

  A /band/ is the open interval of values {(v[k-1] _ v[k])} between
  two consecutive isoline levels. The /index/ of the band is the index {k}
  of its upper level.  */
  
void pswr_compute_band_indices
  ( double vStart, 
    double vStep, 
    double fMin, 
    int kMin, 
    int *iMin, 
    double *zMin, 
    double fMax, 
    int kMax, 
    int *iMax, 
    double *zMax 
  );
  /* Returns the range {iMin .. iMax} of indices of all bands that
    have non-trivial intersection with the value interval 
    {[fMin _ fMax]}. Requires {fMin <= fMax}.  
    
    Considers only isolines with indices in the range {kMin..kMax}, so
    on output we have {kMin <= iMin <= iMax <= kMax+1}; and also
    {v[iMin-1] <= fMin < v[iMin]} and {v[iMax-1] < fMax <= v[iMax]},
    provided that we redefine {v[kMin-1] = -oo} and {v[kMax+1] = +oo}.
    The values {v[iMin]} and {v[iMax]} are returned in {zMin} and {zMax}.
    
    However, when {fMin} and {fMax} coincide with some isoline level
    {v[k]} in {kMin..kMax}, that condition cannot be satisfied, and
    the output is {iMin = iMax = k}  and {zMin = zMax = fMin = fMax} 
    by convention.  */

/* 
  PLOTTING ISOLINES AND BANDS IN A TRIANGLE

  The procedures in this section are given a triangle with corners
  {a=(xa,ya)}, {b=(xb,yb)}, {c=(xc,yc)}, and the respective function
  values {fa}, {fb}, {fc}. Each procedure interpolates an affine
  function through that data, and splits the triangle into bands
  delimited by isolines. */

void pswr_isolines_in_triangle
  ( PSStream *ps,
    double xa, double ya, double fa,
    double xb, double yb, double fb,
    double xc, double yc, double fc,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax        /* Maximum isoline index. */
  );
  /* Plots the isolines that enter the triangle. Specifically, a
    straight isoline segment is drawn, with the current pen, wherever
    the interpolated function matches some isoline level {v[k]}, with
    {k} in {kMin .. kMax}. If {fa = fb = fc = v[k]} for some {k}, then
    all three sides of the triangle are drawn. */
 
void pswr_bands_in_triangle
  ( PSStream *ps,
    double xa, double ya, double fa,
    double xb, double yb, double fb,
    double xc, double yc, double fc,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax,       /* Maximum isoline index. */
    double *R, double *G, double *B
  );
  /* Paints the interior of the triangle with bands whose colors are
    taken from arrays {R[], G[], and B[]}. Specifically, the band
    where the interpolated function lies between isolines {k-1} and
    {k} is painted with color {(R,G,B)[k - kMin]}.
    
    In particular, color {(R,G,B)[0]} is used for the bands where the
    function is less than {v[kMin]}, and color {(R,G,B)[N-1]} for those
    bands where the function is greater than level {v[kMax]}; where 
    {N = kMax - kMin + 2}. */
 
/* 
  PLOTTING ISOLINES AND BANDS IN A QUADRILATRAL

  These procedures are given the corners {(x00,y00), (x01,y01),
  (x01,y10), (x11,y11)} of a quadrilateral, and their respective
  function values {f00,f01,f10,f11}.

  The corners must be given in row-by-row order (NOT ccw order).
  These procedures decompose {Q} into four triangles, each
  determined by one side of {Q} and {Q}'s barycenter {b}. The
  procedures assume that the function's value at {b} is
  {(f00+f01+f10+f11)/4}. The quadrilateral must be star-shaped as
  seen from its barycenter. */

void pswr_isolines_in_quadrilateral
  ( PSStream *ps,
    double x00, double y00, double f00,
    double x01, double y01, double f01,
    double x10, double y10, double f10,
    double x11, double y11, double f11,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax        /* Maximum isoline index. */
  );
  /* Plots isolines in a quadrilateral {Q}.  See {pswr_isolines_in_triangle}
    for more details. */
 
void pswr_bands_in_quadrilateral
  ( PSStream *ps,
    double x00, double y00, double f00,
    double x01, double y01, double f01,
    double x10, double y10, double f10,
    double x11, double y11, double f11,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax,       /* Maximum isoline index. */
    double *R, double *G, double *B
  );
  /* Paints color bands in a quadrilateral {Q}.  See {pswr_bands_in_triangle}
    for more details. */

/* 
  LOW-LEVEL GEOMETRIC TOOLS */

void pswr_sort_triangle
  ( double *xa, double *ya, double *fa,
    double *xb, double *yb, double *fb,
    double *xc, double *yc, double *fc
  );
  /* Permutes the three points {(xa,ya),(xb,yb),(xc,yc)} 
    and the corresponding function values {fa,fb,fc}
    so that {fa <= fb <= fc}. */

void pswr_compute_zero_line_in_triangle
  ( double xa, double ya, double fa,
    double xb, double yb, double fb,
    double xc, double yc, double fc,
    double *xu, double *yu,
    double *xv, double *yv
  );
  /* Given a triangle {a=(xa,ya)}, {b=(xb,yb)}, {c=(xc,yc)}, and
    function values {fa}, {fb}, {fc} at its corners, interpolates an
    affine function through that data, and computes the endpoints
    {u=(xu,yu)} and {v=(xv,yv)} of the line segment within that 
    triangle where the function is zero.
    
    Requires {fa <= fb <= fc} and {fa < fc}. The point {u} will be on the polygonal
    {a--b--c}, with ties broken towards {b}, and the point {v} will be on
    the segment {a--c}. */

#endif
