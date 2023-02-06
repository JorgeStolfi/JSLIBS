#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <epswr.h>
#include <frgb.h>
#include <tf_camera.h>
#include <tf_errors.h>
#include <tf_targets.h>
#include <tf_calib.h>

void tf_plot_cameras
  ( char *out_dir, 
    char *name, 
    int32_t np,
    r3_t p[], 
    int32_t ncpars,
    tf_camera_params_t *cpar[],
    r2_t q[], 
    double hSize,
    double vSize,
    double hMrg,
    double vMrg,
    double Nx,
    double Ny,
    double mRadius,
    int32_t mStyle
  ); 
  /* For each {ip} in {0..ncpars}, if {p} and {cpar[ip]} are not
    NULL, computes the positions {q} from {p} using each
    {cpar[ip]} and plots them with style {mstyle+ip}. Then, if {q}
    is not NULL, plots the marks {q} with style {mStyle + ncpars}.
    All marks are printed with size {mRadius}. */

void tf_plot_marks
  ( epswr_figure_t *eps, 
    int32_t np,
    r2_t q[], 
    double mRadius,
    int32_t mStyle
  );
  /* Plots the image points {q[0..np-1]} to the Postscript stream {eps},
    with style {mStyle}, size {mRadius}, using the current pen 
    and fill color.  The styles are: 0 = circle, 1 = cross,
    2 = diag.cross,  3 = asterisk, 4 = square, 5..oo = diamond. */

