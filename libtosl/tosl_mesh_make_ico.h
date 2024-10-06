/* Creating a keg-shaped {tosl.h} mesh.  */
/* Last edited on 2024-10-06 16:47:42 by stolfi  */

#ifndef tosl_mesh_make_ico_H
#define tosl_mesh_make_ico_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>
#include <tosl_mesh.h>

tosl_mesh_t *tosl_mesh_make_ico(double R);
  /* 
    Creates a mesh data structure for testing the {tosl_mesh_slice}.  
    
    The object is an almost-regural icosahedron with apothem (edge midpoint)
    radius about {R}.  The actual radius may be larger to avoid excessive distortion. */

#endif
