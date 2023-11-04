#ifndef hr2_pmap_opt_translation_H
#define hr2_pmap_opt_translation_H

/* Last edited on 2023-10-20 15:15:27 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>

void hr2_pmap_opt_translation_encode(hr2_pmap_t *M, double rad[], double y[]);
  /* Converts a translation map {M} to two parameters {y[0],y[1]}.
  
    More precisely, assumes that {M} is a translation map with {M.dir[0][0] =
    1}, and that {y}, {rad} have at least 2 elements each. Sets {y[0]}
    and {y[1]} to {M.dir[0][1]/rad[0]} and {M.dir[0][2]/rad[1]},
    respectively. as {y[0]} and {y[1]}of the with {M.dir[0][0] = 1}. */

void hr2_pmap_opt_translation_decode(double y[], double rad[], hr2_pmap_t *M);
  /* Converts two parameters {y[0],y[1]} to a translation map {M}, 
    as the inverse of {hr2_pmap_opt_translation_encode}.
  
    More precisely, assumes that {y}, {rad} have at least 2 elements each,
    and {M} is a translation map with {M.dir[0][0]=M.inv[0][0]=1}. Sets {M.dir[0][1]}
    and {M.dir[0][2]} to {y[0]*rad[0]} and {y[1]*rad[1]}, respectively;
    and sets {M.inv} to be the inverse of {M.dir}. */

#endif
