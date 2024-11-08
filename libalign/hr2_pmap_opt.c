/* See {hr2_pmap_opt.h}. */
/* Last edited on 2024-11-03 16:56:55 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <rn.h>
#include <r3x3.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <affirm.h>
#include <minn.h>
#include <minn_plot.h>
#include <sve_minn.h>

#include <hr2_pmap_encode.h>

#include <hr2_pmap_opt.h>
  
double hr2_pmap_opt_homo_scaling_bias(hr2_pmap_t *M)
  { double sumA2 = r3x3_norm_sqr(&(M->dir));
    double sumB2 = r3x3_norm_sqr(&(M->inv));
    return 0.25*(sumA2 + 1/sumA2 + sumB2 + 1/sumB2);
  }
