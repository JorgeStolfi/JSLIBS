/* The topological slicing algorithm of Minetto et al. (2024).  */
/* Last edited on 2024-10-06 18:53:09 by stolfi  */

#ifndef tosl_mesh_slice_H
#define tosl_mesh_slice_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_slice.h>

/* 
  REFERENCE

    Ricardo Dutra da Silva, Henrique Romaniuk Ramalho, Neri Volpato,
    Rodrigo Minetto, Jorge Stolfi (2024): "An optimal Z-monotone 3D mesh
    slicing algorithm for additive manufacturing". To appear.

  MESH REPRESENTATION

    The algorith in the paper assumes that the mesh to be slices is
    represented by a half-edge data structure, as defined in {haf.h}.
    For efficiency, this implementation assumes a simplified version of
    the half-edge data structure that uses array elements instead of
    heap-allocated nodes, and integer indices instead of pointers.

    Specifically, a mesh with {NE} edges and {NV} vertices is
    represented by a table {Arc[0..2*NE-1} of {tosl_arc_t} records,
    that represent the arcs (oriented edges), and a table
    {Vpos[0..NV-1]} of integer triplets ({tosl_point_t}), which are the
    (quantized) coordinates of the vertices.

    Each {haf_arc_t} {a} of the half-edge is represented by a
    {tosl_arc_t} record {Arc[ia]}, where {ia} is {haf_arc_id(a)}.
    Thus the two arcs of the (unoriented) edge {e} are represented by
    adjacent entries {Arc[2*ke]}} and {Arc[2*ke+1]}, where {ke =
    haf_edge_id(e)}.
  */
      
typedef void tosl_post_proc_t(tosl_plane_id_t ip, tosl_slice_t *S);
  /* Type of a procedure that consumes a slice {S} of a mesh
    by the slicing plane {ip}. */

void tosl_mesh_slice
  ( tosl_mesh_t *mesh,
    int32_t NP,
    tosl_coord_t Zplane[],
    tosl_arc_id_t L[],
    tosl_post_proc_t *post,
    double *time_P
  );

  /* Slices the given {mesh} by a list of planes, using the {Topological-Slicing}
    algorithm of Minetto et al. (2024).
    
    The planes are specified by their (quantized) {Z}-coordinates
    {Zplane[0..NP-1]}, which must be strictly increasing. Every
    {Zplane[ip]} must be distinct from every vertex {Z}-coordinate
    {Vpos[iv].c[2]}.
    
    For every arc {a} of the half-edge, let {ia} be {haf_arc_id(a)}. On
    entry, the field {mesh.Arc[ia].skip} must be set to
    {haf_arc_id(haf_lnext(a))}. These {.skip} fields will be replaced by
    shortcuts during the algorithm.
    
    On entry, also, {mesh.Vpos[E[ia].ivorg]} must be the coordinates of the
    origin vertex of arc {a}.
    
    The procedure calls {post(ip,S)} for the slice {S} obtained
    with {{Zplane[ip]}, for {ip} in {0..NP-1}, in this order. 
    
    If {time_P} is not {NULL}, returns in {*time_P} the
    CPU time used by the topological slicing proper, including
    the allocation of each slice record but not the
    calls to {post}. */

#endif
