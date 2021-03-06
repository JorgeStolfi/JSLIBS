#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <pswr.h>
#include <frgb.h>
#include <tf_camera.h>
#include <tf_errors.h>
#include <tf_targets.h>
#include <tf_calib.h>

void tf_plot_cameras
  ( PSStream *ps, 
    char *name, 
    int np,
    r3_t p[], 
    int ncpars,
    tf_camera_params_t *cpar[],
    r2_t q[], 
    double hSize,
    double vSize,
    double hMrg,
    double vMrg,
    double Nx,
    double Ny,
    double mRadius,
    int mStyle
  ); 
  /* For each {ip} in {0..ncpars}, if {p} and {cpar[ip]} are not
    NULL, computes the positions {q} from {p} using each
    {cpar[ip]} and plots them with style {mstyle+ip}. Then, if {q}
    is not NULL, plots the marks {q} with style {mStyle + ncpars}.
    All marks are printed with size {mRadius}. */

void tf_plot_marks
  ( PSStream *ps, 
    int np,
    r2_t q[], 
    double mRadius,
    int mStyle
  );
  /* Plots the image points {q[0..np-1]} to the Postscript stream {ps},
    with style {mStyle}, size {mRadius}, using the current pen 
    and fill color.  The styles are: 0 = circle, 1 = cross,
    2 = diag.cross,  3 = asterisk, 4 = square, 5..oo = diamond. */

