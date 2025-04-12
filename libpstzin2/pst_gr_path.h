/* Last edited on 2025-03-14 06:50:24 by stolfi */
/* Geometry of a curved line that represents a graph edge in plots. */ 

#ifndef pst_gr_path_H
#define pst_gr_path_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>
#include <r2.h>

typedef struct pst_gr_path_t
  { uint32_t n;     /* number of internal vertices*/
    r2_t* v;        /*coordinates of internal vertices*/
    bool_t reverse; /* TRUE if the vertices are stored in reversed manner*/
  } pst_gr_path_t;
  /* Descriptor of a directed polygonal line on the plane.   The initial and final 
    vertices are defined externally.   The internal vertices are {v[0..n-1]}.
  
    If {reverse} is false, the vertices should be taken in the order {v[0]}
    to {v[n]}.  If {reverse} is true, the vertex sequence {v[0..n-1]} should be
    taken in the reverse order, from {v[n-1]} to {v[0]}. */

#define pst_gr_path_NULL ((pst_gr_path_t){ .n = 0, .v = NULL, .reverse = FALSE })
  /* A null value for a {pst_gr_path_t} variable. */

#define pst_gr_path_MAX_VERTICES (4*4096)
  /* Should be far more than enough for most cases. */

void pst_gr_path_free(pst_gr_path_t P);
  /* Reclaims the storage used by {P} (but not the descriptor {P} itself). */

pst_gr_path_t pst_gr_path_reverse(pst_gr_path_t P);
  /* Returns the path {r} which is the reverse of path {P}.  
    NOTE: {r} and {P} will share the same vertex list (that is, {r.v == P.v}),
    so at most one of them should be reclaimed. */

void pst_gr_path_ctr_dir(r2_t *org, pst_gr_path_t P, r2_t *dst, r2_t *ctr_P, r2_t *dir_P);
  /* Returns the midpoint {ctr} of the path {P} and a vector {dir} parallel to 
    the mean direction of the path around that point.  Returns those values in 
    {*ctr_P} and {*dir_P}.

    More precisely, 
      
      if {p.n == 0, {ctr} is the midpoint of {org} and {dst}, and {dir} is {dst-org};
      
      if {p.n == 1, {ctr} is {P.v[0]} and {dir} is {dst-org};
      
      otherwise, if {P.n} is odd, {ctr} is the middle vertex of {P.v},
      and {dir} is the difference between the nex and previous vertices;
      
      otherwise, {ctr} is the midpoint of the two middle vertices of {P.v},
      and {dir} is their difference. */

r2_t pst_gr_path_start_dir(r2_t *org, pst_gr_path_t P, r2_t *dst);
  /* Returns the initial direction of the path {P}, assuming that 
    it starts at {org} and ends at {dst}.  */

pst_gr_path_t pst_gr_path_concatenate(pst_gr_path_t P0, r2_t *mid, pst_gr_path_t P1);
  /* Returns a path that is the concatenation of {P0} and {P1}, with the
    point {mid} inserted between them.
    
    Typically, {P0} and {P1} describe the plot curves for two arcs, one
    that ends at some vertex with coordinates {mid} and another that
    starts at that vertex. */

void pst_gr_path_write(FILE* wr, pst_gr_path_t P);
  /* Writes to {wr} a decription of the path {P}, in a format compatible with
    {pst_gr_read_path}. */

pst_gr_path_t pst_gr_path_read(FILE* rd);
  /* Reads from {rd} a description of a {pst_gr_path_t}. It should be an integer {n},
    a reversal bit {rev} (0 or 1), and {n} pairs "({X},{Y})" where {X}
    and {Y} are {double} coordinates. */

#endif
