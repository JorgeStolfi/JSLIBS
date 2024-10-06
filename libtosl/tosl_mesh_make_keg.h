/* Creating a keg-shaped {tosl.h} mesh.  */
/* Last edited on 2024-10-06 16:47:53 by stolfi  */

#ifndef tosl_mesh_make_keg_H
#define tosl_mesh_make_keg_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>
#include <tosl_mesh.h>

tosl_mesh_t *tosl_mesh_make_keg
  ( int32_t NS,
    int32_t NR,
    int32_t NB,
    tosl_coord_t Zmax
  );
  /* 
    Creates a mesh data structure for testing the {tosl_mesh_slice}.  
    
    The object is a convext barrel-like solid with {NS}-fold rotational
    symmetry about the {Z}-axis and mirror symmetry about the {Z=0}
    plane. It has horizontal faces at top and botton which are regular
    {NS}-sided polygons. The wall of the barrel is a grid of {NS*NR}
    faces, with consisting of {NR} rings each with {NS} flat and weakly
    convex faces. Each face is a trapezoid with two horizontal sides and
    two ``meridional'' sides. Each of the latter is subdivided into
    {NB+1} collinear edges by {NB} vertices inserted at equal intervals.
    
    The object spans in {Z} a bit beyond  the range {-Zmax .. +Zmax}. */

#endif
