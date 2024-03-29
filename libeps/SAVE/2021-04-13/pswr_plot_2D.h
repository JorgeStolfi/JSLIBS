#ifndef pswr_plot_2D_H
#define pswr_plot_2D_H

/* Plotting of a 2D B�zier patch with fixed resolution (for debugging) */
/* Last edited on 2021-04-13 22:26:03 by jstolfi */

#include <bool.h>
#include <interval.h>

#include <pswr.h>

typedef void pswr_plot_2D_func_t(double x[], int nx, double f[], int nf);
  /* Type of a generic function to be plotted. It should take the
    argument vector {x[0..nx-1]} and store the corresponding function
    value into {f[0..nf-1]}. */
    
typedef struct pswr_plot_2D_style_t
  { int ic; /* Index of first function value to use for smooth color shading. */
    int nc; /* Number of function values to use for smooth color shading. */
    int ib; /* Index of function value to use for discrete color banding. */
    int iv; /* Index of function value to use for isoline drawing. */
    double vStart;
    double vStep;
    int kMin;
    int kMax;
    double *Rtb;
    double *Gtb;
    double *Btb;
  } pswr_plot_2D_style_t;
  /* 
    Type of a record that determines how a two-argument function {F}
    is to be used when filling and/or drawing a shape.
  
    The function {F} is assumed to take a point {x[0..1]} of the plane
    and return a vector {f[0..nf-1]} of values. 
    
    Filling a shape with this style means:

      * if {nc > 0}, {ic >= 0}, and {ic+nc <= nf}, the region is filled
        with smoothly varying colors determined by the function values
        {f[ic..ic+nc-1]}. Each triangular fragment is painted using
        {pswr_shade_triangle} with option {ns == -1} (Gouraud
        shading). The interpretation of the function values depends on
        the number {nc}. If {nc == 1}, the face is painted with a
        monochromatic blue-white-orange scale. If {nc == 2} it is
        painted with a bicolor scale where {(0,0)} is white and the
        two values are orthogonal chroma vectors. If {nc >= 3} the
        first three values are interpreted as RGB values. The
        fields {vStart,vStep,kMin,kMax,Rtb,Gtb,Btb} are ignored.

      * otherwise, if {ib >= 0} and {ib < nf}, the region is filled
        with discrete bands determined by the function value {f[ib]},
        the isoline fields {vStart,vStep,kMin,kMax}, and the color
        tables {Rtb,Gtb,Btb}. See {pswr_bands_in_triangle},
        {pswr_bands_in_quadrilateral}. If these
        tables are NULL, the bands are supressed and the region remains
        unpainted.

      * otherwise, the region is filled with the current fill color.
        The fields {vStart,vStep,kMin,kMax,Rtb,Gtb,Btb}
        are ignored.
        
    Drawing the shape with this style means:
      
      * if {iv > 0} and {iv < nf}, isolines are
      drawn according to the function value {f[iv]} and the isoline
      fields {vStart,vStep,kMin,kMax}. The fields {Rtb,Gtb,Btb}
      are ignored.  See {pswr_isolines_in_triangle},
      {pswr_isolines_in_quadrilateral}.
  */
    
void pswr_plot_2D_tri
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    double xa[],
    double xb[],
    double xc[],
    int ns,
    pswr_plot_2D_style_t *st,
    bool_t fill,
    bool_t draw
  );
  /* Paints a 2D region whose shape and colors are defined by the
    function {F}. The domain of {F} is the triangle {T = (xa,xb,xc)}
    of {R^2}. Each corner of {T} is a vector with two elements.
    
    The domain triangle {T} is divided into {ns^2} triangular
    chips,the function {F} is evaluated at the corners of each chip,
    and each chip is painted as if the function values (shape and
    colors) were affine (1st degree) within the chip.
    
    The function {F} is assumed to take an argument {x[0..1]} and
    return a result vector {F(x) = f[0..nf-1]} with {nf} elements. The
    first 2 components of {f[0..nf-1]} are interpreted as the plotting
    coordinates {P(x)} of the point {x}. The remaining coordinates
    (i.e., {f[2..nf-1]} are displayed by means of colors and/or
    isolines.
    
    If {fill} is true, the region is first filled as specified by the
    {st} and the current fill color.
    
    After filling (if any), if {draw} is TRUE, the region is overlaid
    with lines drawn as specified by the {st} and the current pen
    settings. */

void pswr_plot_2D_tri_outline
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    double xa[],
    double xb[],
    double xc[],
    int ns
  );
  /* Draws (with the current pen settings) the outline of a 2D region
    whose shape is defined by the function {F} and the triangle {T = (xa,xb,xc)},
    as explained for  {pswr_plot_2D_tri}.
    
    The parameters have the same meaning as in {pswr_plot_2D_tri}.
    Each side of the domain triangle {T} is divided into {ns}
    segments, the endpoints of those segments are mapped through {F}
    to points of the plane, and these points are connected by straight
    lines. */

void pswr_plot_2D_quad
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    interval_t B[],
    int n0,
    int n1,
    pswr_plot_2D_style_t *st,
    bool_t fill,
    bool_t draw
  );
  /* Paints a 2D region whose region and colors are defined by the
    function {F}. The domain of {F} is the the rectangle 
    {B = B[0]�B[1]} of {R^2}.
    
    The domain rectangle {B} is divided into {n0 � n1} rectangular
    tiles, each tile is divided into 4 triangular chips, the function
    {F} is evaluated at the corners of each chip, and each chip is
    painted as if the function values (shape and colors) were affine
    (1st degree) within the chip.
    
    The other parameters have the same meaning as those of 
    {pswr_plot_2D_tri}. */

void pswr_plot_2D_quad_outline
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    interval_t B[],
    int n0,
    int n1
  );
  /* Draws (with the current pen settings) the outline of a 2D region
    whose shape is defined by the function {F} and the 
    rectangle {B = B[0]�B[1]}, as explained for  {pswr_plot_2D_quad}.
    
    The parameters have the same meaning as in {pswr_plot_2D_quad}.
    Each side of the domain rectangle {B} is divided into {n0} or {n1}
    segments, according to the corresponding axis of the domain; the
    endpoints of those segments are mapped through {F} to points of
    the plane, and these points are connected by straight lines. */

void pswr_plot_2D_line
  ( PSStream *ps, 
    int nf,
    pswr_plot_2D_func_t *F,
    interval_t *B,
    int ns
  );
  /* Draws (with the current pen settings) a curved line 
    whose shape is defined by the function {F}, that maps    
    the real interval {*B} of {R} into {R^nf}.
    
    The domain interval {*B} is divided into {ns} sub-intervals;
    the function is evaluated at the ends of each sub-interval, 
    and the resulting points are connected by straight lines,
    drawn with the current pen settings.
    
    The function {F} is assumed to take a single-element vector
    argument {x[0..0]} and return a result vector {F(x) = f[0..nf-1]}
    with {nf} elements. The first 2 components of {f[0..nf-1]} are
    interpreted as the plotting coordinates {P(x)} of the point {x}.
    The remaining coordinates (i.e., {f[2..nf-1]} are ignored. */

/* 
  AFFINE PAINTING PROCEDURES

  These procedures plot the given triangle as if the function values
  were bilinear inside it. 
  
  The function values at the corners {fa,fb,fc} must be given. Each
  argument must be a vector with {nf} values. Function values {f[0]}
  and {f[1]} are used as the X and Y coordinates for plotting. */
  
void pswr_plot_2D_tri_atom
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf,
    pswr_plot_2D_style_t *st,
    bool_t fill,
    bool_t draw
  );
  /* If {fill} is true, fills the interior of the triangle with
    colors; then, if {draw} is true, draws the isolines with
    the current pen. The color interpretation and the remaining
    parameters have the same meanings as in {pswr_plot_2D_quad}. */

void pswr_plot_2D_tri_atom_solid
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf 
  );
  /* Paints the triangle with the current fill color. */

void pswr_plot_2D_tri_atom_shade
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf,
    int ic,
    int nc
  );
  /* Paints the triangle with smooth shading. The colors at the
    vertices are defined by the function values {f[ic..ic+nc-1]}. */
    
void pswr_plot_2D_tri_atom_bands
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf,
    int ib,
    double vStart,
    double vStep, 
    int kMin,
    int kMax,
    double *Rtb, 
    double *Gtb, 
    double *Btb 
  );
  /* Paints the triangle {fa,fb,fc} with discrete color bands defined
    by the the function value {f[ib]}, the isoline parameters
    {vStart,vStep,kMin,kMax}, and the color tables {Rtb,Gtb,Btb}. If
    the tables are NULL, paints nothing. If they are non-NULL, each
    must have {N = kMax - kMin + 2} elements. Assumes that {f[ib]}
    varies affinely inside the triangle. See
    {pswr_bands_in_triangle}. */

void pswr_plot_2D_tri_atom_isolines
  ( PSStream *ps, 
    double fa[],
    double fb[],
    double fc[],
    int nf,
    int iv,
    double vStart,
    double vStep, 
    int kMin,
    int kMax 
  );
  /* Draws isolines inside the triangle {fa,fb,fc}
    with the current pen settings.
    The isolines are defined by the function value {f[iv]} and the
    parameters {vStart,vStep,kMin,kMax}.  Assumes that {f[ib]} varies
    affinely inside the triangle. See
    {pswr_isolines_in_triangle}. */

#endif
