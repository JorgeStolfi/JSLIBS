#ifndef hxg_eps_H
#define hxg_eps_H

/* Routines to plot hexagonal canvas pixels and other things. */
/* Last edited on 2025-01-01 01:45:33 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <r2.h>
#include <interval.h>
#include <i2.h>
#include <epswr.h>
#include <frgb.h>

#include <hxg_canvas.h>

epswr_figure_t *hxg_eps_new_figure(interval_t B[], char *fname);
  /* Creates a figure object that writes to file {fname}, and sets up
    the scales and windows for plotting things in the ranges
    {B[0],B[1]}. Use the commands in {epswr.h} or {hxg_eps.h} for
    plotting. Be sure to call {epswr_end_figure} before exiting the
    program. */

void hxg_eps_plot_canvas
  ( epswr_figure_t *eps,  /* Where to plot. */
    hxg_canvas_t *cvs,    /* The canvas to plot. */
    double pixv[],        /* Pixel values. */
    uint32_t maxv,        /* Max pixel value in tables. */
    double *r,            /* Radius table. */
    frgb_t *fill_color,           /* Fill color table. */
    frgb_t *draw_color            /* Draw color table. */
  );
  /* Plots the pixels of the canvas {cvs} into the Postscript file
    {eps}.  
    
    The array {pixv} must have {cvs.nx*cvs.ny} elements, and has the
    grid pixel values stored by rows. The arraus {r,fill_color,draw_color} must have
    {maxv+1} elements each. A pixel whose value rounds to an integer {v}
    in the range {0..maxv} is plotted as a circle with radius {r[v]}
    filled with color {fill_color[v]} and drawn with color {draw_color[v]}.
    
    If the rounded pixel value {v} is negative, the absolute value of
    {v} is used, but the fill color is complemented. If the rounded value
    {v} exceeds {maxv}, then {maxv} is used. */

#endif
