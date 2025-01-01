#ifndef cpk_eps_H
#define cpk_eps_H

/* Plot functions for {libcpk}. */
/* Last edited on 2024-12-31 20:21:22 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <r2.h>
#include <interval.h>
#include <i2.h>
#include <epswr.h>
#include <frgb.h>

epswr_figure_t *cpk_eps_new_figure(interval_t B[], char *fname);
  /* Creates a figure object that writes to file {fname}, and sets up
    the scales and windows for plotting things in the ranges
    {B[0],B[1]}. Use the commands in {epswr.h} or {cpk_eps.h} for
    plotting. Be sure to call {epswr_end_figure} before exiting the
    program. */

void cpk_eps_fill_polygon(epswr_figure_t *eps, r2_vec_t *p);
  /* Fills the polygon {p} using the current fill color. */

void cpk_eps_draw_polyline(epswr_figure_t *eps, r2_vec_t *p);
  /* Draws the polyline(s) whose vertices are {p} 
    using the current pen width and dash settings, and the
    current draw color.  */

void cpk_eps_plot_circle
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_t *c,             /* Center of circle. */
    double r,            /* Radius of circle. */
    bool_t fill, bool_t draw /* Plotting options. */
  );
  /* Fills and/or draws the circle with center {c} and radius
    {r} using the current pen width and dash settings, and the
    current fill and draw colors. */

void cpk_eps_plot_dot
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_t *c,             /* Center of circle. */
    double r,            /* Radius of circle (mm). */
    bool_t fill, bool_t draw /* Plotting options. */
  );
  /* Same as {cpk_eps_plot_circle}, except that the 
    radius is in actual millimeters, irrespective of 
    the current scale. */

void cpk_eps_plot_rect_stroke
  ( epswr_figure_t *eps, /* Where to plot. */
    r2_t *a,             /* Start of stroke. */
    r2_t *b,             /* End of stroke. */
    double r,            /* Half-width of stroke. */
    bool_t fill, bool_t draw /* Plotting options. */
  );
  /* Fills and/or draws the sausage that is the circles with centers {a} and {b}
    plus a rectangle connecting the two points, all with radius,
    {r} using the current pen width and dash settings, and the
    current fill and draw colors. */

void cpk_eps_fill_circles
  ( epswr_figure_t *eps, /* Where to plot. */
    r2_vec_t *c,         /* Centers of circles. */
    double_vec_t *h,     /* Individual radii. */
    double r             /* Fixed radius. */
  );
  /* Calls {cpk_eps_plot_circle(eps, c[k], h[k]+r, TRUE, FALSE)}
    for all the {k}.  If {h} is NULL or shorter than {c}, 
    the missing entries are assumed to be 0. */

void cpk_eps_fill_dots
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_vec_t *c,          /* Centers of circles. */
    double r             /* Dot radius (mm). */
  );
  /* Calls {cpk_eps_plot_dot(eps, c[k], r, TRUE, FALSE)}
    for all the {k}. */

void cpk_eps_show_labels
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_vec_t *c,         /* Reference points. */
    r2_t disp,           /* Extra displacement for all ref pts. */
    double xalign,       /* Horiz rel pos of label to align with ref X. */
    double yalign,       /* Vert rel pos of label to align with ref Y. */
    double ptsize        /* Text point size. */
  );
  /* For each {k} in {0..c.ne-1}, writes the string "{k}" 
    at the point {c[k] + disp}, with font of size {ptsize}.
    The parameters {xalign} and {yalign} behave as in {epswr_label}. */
    
void cpk_eps_draw_circles
  ( epswr_figure_t *eps,      /* Where to plot. */
    r2_vec_t *c,        /* Centers of circles. */
    double_vec_t *h,    /* Individual radii. */
    double r           /* Fixed radius. */
  );
  /* Calls {cpk_eps_plot_circle(eps, c[k], h[k]+r, FALSE, TRUE)}
    for all the {k}.  If {h} is NULL or shorter than {c}, 
    the missing entries are assumed to be 0. */

void cpk_eps_fill_sausage
  ( epswr_figure_t *eps,        /* Where to plot. */
    r2_vec_t *p,          /* Vertices of sausage. */
    double r             /* Half-width of sausage. */
  );
  /* Fills the sausage-shape with half-width {r} whose skeleton is the
    polyline {p}. Said another way, strokes the polyline {p} with a
    round pen of radius {r} --- but with the current *fill* color.
    
    Equivalent to {cpk_eps_plot_circle(eps, p[k], r, TRUE, FALSE)} for
    every finite point {p[k]}, and {cpk_eps_plot_rect_stroke(eps, p[k-1],
    p[k], r, TRUE, FALSE)} for all consecutive finite point pairs
    {p[k-1], p[k]}. */

#endif
