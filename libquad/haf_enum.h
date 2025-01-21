#ifndef haf_enum_H
#define haf_enum_H

/* Enumerating and renumbering elements of a half-edge structure. */
/* Last edited on 2025-01-09 23:19:38 by stolfi */

#define half_enum_H_copyright \
  "Copyright (C) 2023 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

#include <haf.h>

void haf_enum_edges
  ( haf_arc_count_t NR,
    haf_arc_t root[],
    haf_edge_id_t eid0,
    haf_edge_vec_t *E_P,
    haf_edge_id_vec_t *oid_P,
    haf_arc_vec_t *C_P
  );
  /* Given a list {root[0..NR-1]} of arc references, enumerates all arc
    references that are reachable from them through any sequence of
    {.sym} and {.lnext} operators.
    
    As a side effect, the procedue reassigns the ids {.eid} of all
    reached edges as sequential indices starting from {eid0}. The two
    opposite arcs of each edge are assigned arc ids {2*eid} and
    {2*eid+1}, in unspecified order.
    
    The procedure builds a table {E} of {haf_edge_t}s and a table {oid}
    of {haf_edge_id_t}s, such that {E.ne} and {oid.ne} will be the the
    number {NE} of edges reached and renumbered; {E.e[ke]}, for {ke} in
    {0..NE-1}, is the edge whose id is {eid0+ke}; and {oid.e[ke]}] will
    be the original ID of the edge {E.d[ke]}, before it was renumbered.
    
    The table will also build a table {C} of arcs, such that {C.NE} will
    the number {NC} of connected components formed by those edges, and
    {C.e[0..NC-1]} is one base arc on each connected component. These
    arcs will be a subset of {root[0..NR-1]} or their {.sym}.
    
    If {E_P} is not {NULL}, the table {E} will be returned in {*E_P}. The
    corresponding array {E_P->e} will be reused, overwritting any previous
    contents, and is resized as needed.  Ditto for {C_P} and {C}, 
    and for {oid_P} and {oid}. */

typedef uint64_t haf_face_id_t; /* Numeric identifier of a face of the structure. */

typedef uint64_t haf_vert_id_t; /* Numeric identifier of a vertex of the structure. */

void haf_enum_faces(haf_edge_count_t NE, haf_arc_t a[], haf_edge_id_t eid0, haf_face_id_t fid[]);
void haf_enum_verts(haf_edge_count_t NE, haf_arc_t a[], haf_edge_id_t eid0, haf_vert_id_t vid[]);
  /* These procedures assume that {a[0..NE-1]} are the base arcs of all
    edges of the. The id numbers of those edges must satisfy 
    {a[ke].eid=ke+eid0} for all {ke} in {0..NE-1}.
    
    The vectors {fid} and {vid} must have {2*NE} elements. For each arc
    {e} in the structure, {haf_get_faces} stores into
    {fid[e.aid-2*eid0]} a unique numeric identifier of the left face
    ({.lnext} cycle) of arc {e}. These face ids are assigned
    sequentially from 0, in the order the faces are found.
    
    The effect of {haf_get_verts} is the same, except that it identifies
    vertices ({.oprev} cycles) among the arcs {a[0..NE-1]} and their
    opposites.
    
    The set of the arcs {a[..NE-1]} and their opposites must be closed
    under the arc walking operators, meaning that the {.lnext} of any
    arc {a[0..NE-1]} and of its opposite must be some arc {a[0..NE-1]}
    or its opposite. Thus the vector {a} can be, for example, {E.e}
    where {E} is the result of {haf_enum_edges(NR, root, eid0)}.*/

#endif
