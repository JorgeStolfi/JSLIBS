#ifndef gem_bary_H
#define gem_bary_H
/* Last edited on 2015-12-01 15:47:35 by stolfilocal */

/* 
  BARYCENTRIC GEMS
  
  A barycentruc gem is the barycentric subdivision
  of some {d}-dimensional map.  In each simplex of the 
  gem, the corner with color {k} lies in 
  some part of the map with dimension {k},
  for {k} in {0,1,2...d}.
  
*/

#define _GNU_SOURCE
#include <gem.h>

void gem_bary_splice(gem_ref_t a, gem_ref_t b, int k);
  /* Glues or splits two {d}-dimensional parts of 
    a {d}-dimensional map by a shared {k}-part (facet),
    where {d=k+1}. 
    
    Either the facets specified by {a} and {b}
    are distinct but isomorphic, and both on the border of the {d}-map, 
    or they are the same interior facet of the {d}-map.  In the first case,
    the two facets are glued and become the same 
    interior facet.  In the second case, the two {d}-parts
    are separated and that shared interior facet
    becomes two border facets.
    
    Assumes that the
    barycentric subdivision of each facet consists of walls with
    wall-color {d+1}, belonging to cells that share one vertex with
    color {d} and are connected through walls with colors {0..d-1}. The
    nodes {a} and {b} must be corresponding cells incident to that face.
    All the nodes must have been allocated with dimension parameter
    {d+1} or more. */

#endif
