#ifndef hr2_pmap_throw_type_H
#define hr2_pmap_throw_type_H

/* Generate random projective maps of a given type and handedness. */
/* Last edited on 2024-09-17 16:29:52 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <hr2_pmap.h>
#include <hr2.h>

hr2_pmap_t hr2_pmap_throw_type(hr2_pmap_type_t type, sign_t sgn);
  /* Returns a random projective map with given {type} and 
    handedness {sgn}. */

#endif
