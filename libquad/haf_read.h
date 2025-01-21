#ifndef haf_read_H
#define haf_read_H

/* Reading the half-edge representation of 2D meshes. */ 
/* Last edited on 2025-01-09 23:42:00 by stolfi */

#define half_read_H_copyright \
  "Copyright (C) 2023 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <stdio.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

#include <haf.h>

haf_arc_t haf_read_arc(FILE *rd, haf_edge_id_t eid0, haf_edge_vec_t *A);
  /* Parses from file {rd} a description of an arc, in the format
    "{eid}:{db}" where {eid} is a numeric edge id and {db} is a
    direction bit. Requires a table {A} such that {A.e[ke]} is the edge
    with id {eid0 + ke}. The resulting arc {a} has {haf_edge_id(a) ==
    {eid}} and {haf_dir_bit(a) == db}. */

void haf_read_map(FILE *rd, haf_edge_vec_t *A, haf_edge_id_t *eid0_P, haf_arc_vec_t *R);
  /* Reads from {rd} a description of a map in the format used by
    {write_map}, and builds the corresponding half-edge structure in
    memory.
    
    The procedure allocates all the edge records specified in {rd}, and
    assigns to them sequential edge id numbers starting from the
    {edge_id_min} parameter found in the file's header (which is
    returned in {*eid0_P}). Then it stores those edges into {A.e[0..A.ne-1]}
    It also stores into {R.e[0..R.ne-1]} the
    root arcs specified in the file. 
    
    The parameters {A} and {R} must not be null and the vectors
    must be initialized by the user.  The arrays {A.e} and
    {R.e} are resized as needed. If they are not empty on entry, 
    their storage reused, and their previous contents is overwritten. */

#endif
