/* Test tools for {r2_opt.h}. */
/* Last edited on 2017-06-06 20:54:05 by stolfilocal */

#ifndef test_r2_opt_tools_H
#define test_r2_opt_tools_H

#define _GNU_SOURCE
#include <stdint.h>
#include <r2.h>
#include <bool.h>

#include <test_r2_opt_basic.h>

double tr2o_tools_eval_mother_image
  ( r2_t *p, 
    r2_t *scale,
    int mom_NF, 
    r2_t mom_frq[], 
    r2_t mom_phi[], 
    double mom_amp[]
  );
  /* An image which is the sum of {mom_NF} sine-grid waves with various
    wavelengths.
    
    Each wave {k} is the product of two sine waves with 
    same absolute frequency but perpendicular normal directions.
    The spatial frequency of the first factor is {mom_frq[k]},
    that of the secod factor is rotated 90 degrees from that.
    The initial phase of each factor {r} (0 or 1) is {mom_phi[k].c[r]}.
    The amplitude of the product is {mom_amp[k]}. Assumes that the tables
    {mom_amp,mom_phi,mom_frq} have been initialized.
    
    The image is internally smoothed so that it has no components with
    wavelength greater than or equal to {2*scale.c[j]} along axis {j}.
    In other words, it is evaluated as if it was shrunk by
    {1/scale.c[j]} along each axis {j}, with proper antialiasing, and
    the point {p} was reduced by the same amount. */

void tr2o_tools_initialize_mother_image
  ( int mom_NF, 
    r2_t mom_frq[], 
    r2_t mom_phi[], 
    double mom_amp[],
    bool_t verbose
  );
  /* Initalizes the tables {mom_amp,mom_phi,mom_frq} for
    {tr2o_tools_mother_image}, each with {mom_NF} elements. */

#endif
