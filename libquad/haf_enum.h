#ifndef haf_enum_H
#define haf_enum_H

/* Enumerating and renumbering elements of a half-edge structure. */
/* Last edited on 2024-12-05 10:38:53 by stolfi */

#define half_enum_H_copyright \
  "Copyright (C) 2023 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

#include <haf.h>

void haf_enum_edges
  ( haf_arc_count_t nr,
    haf_arc_t root[],
    haf_edge_id_t eid0,
    haf_arc_vec_t *E,
    haf_arc_vec_t *C
  );
  /* Given a list {root[0..nr-1]} of arc references, enumerates all arc
    references that are reachable from them through any sequence of
    {.sym} and {.lnext} operators.
    
    As a side effect, reassigns all edge ids {.eid} as sequential
    indices starting from {eid0}. The two opposite arcs of each edge are
    assigned arc ids {2*eid} and {2*eid+1}, in unspecified order.
    
    On return, if {E} is not {NULL}, {ne=E.ne} will be the the number of
    edges reached and renumbered, and {E.e[ke]}, for {ke} in {0..ne-1},
    is the base arc on the edge whose id is {eid0+ke}.
    
    On return, if {C} is not {NULL}, {nc=C.ne} will the number of
    connected components formed by those edges, and {C.e[0..nc-1]} is
    one base arc on each connected component. These arcs will be a subset of
    {root[0..nr-1]} or their {.sym}.  
    
    The arrays {E.e} and {C.e} are resized as needed. If either is not {NULL}
    on entry, it is reused, and their previous contents is overwritten. */

typedef uint64_t haf_face_id_t; /* Numeric identifier of a face of the structure. */

typedef uint64_t haf_vert_id_t; /* Numeric identifier of a vertex of the structure. */

void haf_enum_faces(haf_edge_count_t ne, haf_arc_t a[], haf_edge_id_t eid0, haf_face_id_t fid[]);
void haf_enum_verts(haf_edge_count_t ne, haf_arc_t a[], haf_edge_id_t eid0, haf_vert_id_t vid[]);
  /* These procedures assume that {a[0..ne-1]} are the base arcs of all
    edges of the. The id numbers of those edges must satisfy 
    {a[ke].eid=ke+eid0} for all {ke} in {0..ne-1}.
    
    The vectors {fid} and {vid} must have {2*ne} elements. For each arc
    {e} in the structure, {haf_get_faces} stores into
    {fid[e.aid-2*eid0]} a unique numeric identifier of the left face
    ({.lnext} cycle) of arc {e}. These face ids are assigned
    sequentially from 0, in the order the faces are found.
    
    The effect of {haf_get_verts} is the same, except that it identifies
    vertices ({.oprev} cycles) among the arcs {a[0..ne-1]} and their
    opposites.
    
    The set of the arcs {a[..ne-1]} and their opposites must be closed
    under the arc walking operators, meaning that the {.lnext} of any
    arc {a[0..ne-1]} and of its opposite must be some arc {a[0..ne-1]}
    or its opposite. Thus the vector {a} can be, for example, {E.e}
    where {E} is the result of {haf_enum_edges(nr, root, eid0)}.*/

#endif
