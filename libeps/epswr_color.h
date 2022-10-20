/* Useful color operations for Postscript plotting. */
/* Last edited on 2022-10-20 07:38:22 by stolfi */

#ifndef epswr_color_H
#define epswr_color_H

#include <stdint.h>

void epswr_make_color_table
  ( double vStart,
    double vStep,
    int32_t kMin,
    int32_t kMax, 
    double RMin, double GMin, double BMin,
    double RZer, double GZer, double BZer,
    double RMax, double GMax, double BMax,
    int32_t *NP, 
    double **RP, 
    double **GP,
    double **BP
  );
  /* Returns in {*RP,*GP,*BP} a color table adequate for color-band 
    plotting of bivariate functions. The three vectors are allocated
    by the routine, and the number of entries {N} is returned in {*NP}.
    
    Color {R[i],G[i],B[i]} is meant to be used for the band between
    isolines {kMin+i-1} and {kMin-i}. In particular, color
    {R[0],G[0],B[0]} will be {RMin,GMin,BMin}, and is intended for the
    band below isoline {kMin}. Color {R[N-1],G[N-1],B[N-1]} will be
    {RMax,GMax,BMax}, and is intended for the band above isoline
    {kMax}. Therefore, {N = kMax - kMin + 2}.  The band that brackets
    the value 0, if any, will get color {RZer,GZer,BZer}.
    
    More precisely, the color {R[i],G[i],B[i]} assigned to band {k =
    kMin + i} depends on the mean value {m[k] = (v[k-1]+v[k])/2} of
    that band, where {v[k] = epswr_level(vStart, vStep, k)}. If {m[k]}
    is positive, the color is interpolated between {RZer,GZer,BZer}
    and {RMax,GMax,BMax} with ratio {r = +m[k]/m[kMax+1]}. If {m[k]}
    is positive, the color is interpolated between {RZer,GZer,BZer}
    and {RMin,GMin,BMin} with ratio {r = -m[k]/m[kMin]}. */

void epswr_interpolate_colors
  ( double r, 
    double R0, double G0, double B0,
    double R1, double G1, double B1,
    double *R, double *G, double *B
  );
  /* Interpolates between the colors {R0,G0,B0} and {R1,G1,B1} in the
    ratio {r}. The interpolated color lies on a curved path between
    the two colors, in such a way that equal increments in {r} give
    colors that are approximately equidistant in perceptual
    distance.  The resulting color is clipped to the unit RGB cube. */
 
void epswr_color_scale_1
  ( double fs, 
    double Rs, double Gs, double Bs,
    double Y0,
    double *R, double *G, double *B
  );
  /* Maps a function value {fs} from {[-1 _ +1]} to a color {R,G,B},
    so that {0} maps to gray with intensity {Y0}, {1} maps to {Rs,Gs,Bs},
    and {-1} to the complementary hue with same brightness. */
  
void epswr_color_scale_2
  ( double fs, double ft, 
    double Rs, double Gs, double Bs,
    double Y0,
    double *R, double *G, double *B
  );
  /* Maps two function values {(fs,ft)} to a color {R,G,B}, so that
    {(0,0)} maps to gray with intensity {Y0}, {(1,0)} maps to
    {(Rs,Gs,Bs)}, {(0,1)} maps to an orthogonal color {(R1,G1,B1)},
    and {(-1,0)}, {(0,-1)} map to the complementary colors, all with
    the same brightness.
    
    !!! Improve the formulas so that the valid ranges for {fs,ft} are predictable. !!! */
    
void epswr_color_scale_3
  ( double fs, double ft, double fu, 
    double Ymin,
    double Ymax,
    double *R, double *G, double *B
  );
  /* Maps three function values {fs,ft,fu} in {[0 _ 1]} to a color
    {R,G,B}, so that {(-1,-1,-1)} maps to gray with intensity {Y0},
    and {(+1,+1,+1)} maps to a gray with intensity {Y1}. */
    
#endif
