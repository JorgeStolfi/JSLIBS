#ifndef haf_copy_H
#define haf_copy_H

/* Enumerating and renumbering elements of a half-edge structure. */
/* Last edited on 2025-01-09 23:16:19 by stolfi */

#define half_enum_H_copyright \
  "Copyright (C) 2023 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

#include <haf.h>

void haf_copy
  ( haf_arc_count_t NR,
    haf_arc_t root[],
    haf_edge_id_t eid0,
    haf_edge_vec_t *new_P,
    haf_edge_vec_t *old_P
  );
  /* Given a list {root[0..NR-1]} of arc references, enumerates all arc
    references that are reachable from them through any sequence of
    {.sym} and {.lnext} operators, and creates new copies of their
    edge records, linked in isomorphic ways.
    
    The IDs of the new edges will be assigned sequentially as they are
    found, starting with {eid0}  The two opposite arcs of each edge are
    assigned arc ids {2*eid} and {2*eid+1}, in unspecified order.  The IDs of the 
    original edges are not modified.
    
    The procedure builds a table {new} of {haf_edge_t}s such that
    {new.NE} is the number {NE} of edges reached and copied, and
    {new.e[ke]}, for {ke} in {0..NE-1}, is the new edge whose
    (new) id is {eid0+ke}.
    
    The procedure also builds a table {old} of {haf_edge_t}s, also with
    {NE} elements, such that {old.e[ke]} is the original edge whose copy
    is {new[ke]}.
    
    If {new_P} is not {NULL}, the table {new} is returned in {*new_P}. The
    corresponding array {new_P->e} will be reused, overwritting any previous
    contents, and is resized as needed. Ditto for {old_P} and {old}. */

#endif
