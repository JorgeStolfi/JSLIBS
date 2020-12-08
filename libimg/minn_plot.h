#ifndef minn_plot_H
#define minn_plot_H

/* Plots a minimizer's goal function in a 2D subspace of its domain. */
/* Last edited on 2017-03-13 21:33:57 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>

typedef double minn_plot_goal_t(int n, double x[]);

void minn_plot_2D_float_image
  ( FILE *wr, 
    int n,
    minn_plot_goal_t *F,
    double xo[],
    double xa[],
    double xb[],
    int NS
  );
  /* Plots the value of {F(n,x)} for several vectors {x[0..n-1]}.
    Each vector is obtained as {xo + a*(xa-xo) + b*(xb-xo)},
    where ech parameter {a,b} assumes {2*NS+1} equally spaced 
    values in {[-1_+1]}. 
    
    The output is written to {wr} as an image in the
    FNI format, with 1 channel and {2*NS+1} rows and columns. */

void minn_plot_2D_gnuplot(FILE *wr, int n, minn_plot_goal_t *F, double x0[], int NS, double R);
  /* Writes to {wr} a file with sampled values of
    {F} for various parameter vectors {x[0..n-1]}, in a format
    suitable for the {splot} command of {gnuplot}.
    
    If {nx == 1}, the parameter vectors are obtained by setting
    {x[0]} to {2*N+1} equally spaced values spanning the
    interval{x0[0] + [-R_+R]}. If {nx >= 2), the sample vectors {x}
    form a square two-dimensional grid in {R^nx} with {2*N+1} nodes
    and side {2*R} on each axis, centered at the configuration
    {x0[0..nx-1]}. If {nx==2}, the two grid axes are {x[0]} and
    {x[1]}; otherwise they are two random orthogonal directions in
    {R^{nx}}. */
    
void minn_plot_1D_gnuplot(FILE *wr, int n, minn_plot_goal_t *F, double x0[], int NS, double r);
  /* Writes to {wr} the values {F(n,x0+e*u)} for various directions
    {u} in {\RR^{n}} and for {2*NS+1} equally spaced values of {e}
    spanning {[-r _ +r]}.  The directions {u} include the cardinal 
    directions. */

#endif
