/* Translation maps in 2D oriented projective plane {\RT^2}. */
/* Last edited on 2024-12-05 10:27:13 by stolfi */ 
   
#ifndef hr2_pmap_translation_H
#define hr2_pmap_translation_H

#include <stdint.h>

#include <r2.h>
#include <hr2.h>
#include <hr2_pmap.h>

hr2_pmap_t hr2_pmap_translation_from_disp(r2_t *disp);
  /* Returns a projective map that performs a translation by the
    Cartesian vector {disp}. The map is a special case of affine
    (has {[1 0 0]} as the first column), similarity, and congruence. */

#endif
