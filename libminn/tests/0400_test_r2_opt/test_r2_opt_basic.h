/* Basic test defs for {r2_opt.h}. */
/* Last edited on 2025-03-19 14:36:11 by stolfi */

#ifndef test_r2_opt_basic_H
#define test_r2_opt_basic_H

#include <stdio.h>
#include <stdint.h>

#include <r2.h>
#include <i2.h>
#include <bool.h> 

typedef double tr2o_image_eval_proc_t(uint32_t i, i2_t iscale, double x, double y); 
  /* Type of a procedure that evaluates image number {i} at the
    domain point {(x,y)}. The image's domain is implicitly shrunk by
    {1/2^iscale.c[j]} along each axis {j} (with antialiasing).  The point {(x,y)} is
    not scaled. */

void tr2o_debug_params
  ( uint32_t indent,
    char *pname,
    uint32_t NI, 
    r2_t p[]
  );
  /* Prints a list of {r2_t} valued parameters {p[0..NI-1]}.  Indents by {indent}
    columns. */

void tr2o_debug_r2
  ( uint32_t indent,
    char *pname,
    uint32_t i,
    r2_t *p
  );
  /* Prints the point {p} with label "{pname}[{i}]". Indents by {indent} columns. */

void tr2o_debug_points
  ( uint32_t indent, 
    char *title, 
    uint32_t NI, 
    char *pname, 
    r2_t p[], 
    r2_t pini[], 
    r2_t popt[], 
    r2_t arad[],
    r2_t astp[],
    double f2p, 
    double b2p
  );
  /* Prints the points {p[0..NI-1]}, named {pname},
    and the corresponding raw goal function value {f2p}
    and bias term {b2p}.  
    
    If {pini} is not NULL, prints also the difference between 
    each {p[i]} and the corresponding {pini[i]}.  Ditto if {popt} is not NULL.
    The differences are expressed also as multiples of {arad} and {astp}.
    Indents by {indent} columns. */

#endif
