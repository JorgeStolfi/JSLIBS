#ifndef hr2_pmap_opt_H
#define hr2_pmap_opt_H

/* Tools for optimizing projective maps. */
/* Last edited on 2024-11-03 16:57:08 by stolfi */ 

#define _GNU_SOURCE

#include <bool.h>
#include <stdint.h>
#include <r3x3.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <i2.h>

typedef double hr2_pmap_opt_func_t(hr2_pmap_t *M); 
  /* Type of a procedure that evaluates some badness function of the
    projective map {*M}. The function had better achieve a minimum value
    for some finite {*M}. */

typedef bool_t hr2_pmap_opt_pred_t(hr2_pmap_t *A, double fA);
  /* Type of a procedure that detgermines whether the map {A},
    whose badness is {fA}, is acceptable for some purpose. */

double hr2_pmap_opt_homo_scaling_bias(hr2_pmap_t *M);
  /* Returns {f2(M.dir)+f2(M.inv)}, where {f2} is a function of a 3x3
    matrix that is minimum when the sum of the squares of the elements
    is 1.
    
    This function can be used as a bias term when optimizing a general
    projective map, to remove the spurious degree of freedom due to
    homogeneous scaling of the matrices {M.dir} and {M.inv}. */
 
#endif
