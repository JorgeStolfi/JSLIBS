#ifndef haf_read_H
#define haf_read_H

/* Reading the half-edge representation of 2D meshes. */ 
/* Last edited on 2023-10-05 20:25:46 by stolfi */

#define half_read_H_copyright \
  "Copyright (C) 2023 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>

#include <haf.h>

haf_arc_t haf_read_arc(FILE *rd, haf_edge_id_t eid0, haf_arc_vec_t *A);
  /* Parses from file {rd} a description of an arc, in the format
    "{eid}:{db}" where {eid} is a numeric edge id and {db} is 
    a direction bit. Requires a table {A}, which is a vector of arcs such
    that {A.e[ke]} is the base arc on the edge {ed} with edge it {eid0 + ke}.
    The resulting arc {a} has {haf_edge_id(a) == {eid}} and
    {haf_dir_bit(a) == db}. */

void haf_read_map(FILE *rd, haf_arc_vec_t *A, haf_edge_id_t *eid0_P, haf_arc_vec_t *R);
  /* Reads from {rd} a description of a map in the format used by
    {write_map}, and builds the corresponding quad-edge structure in
    memory.
    
    The procedure allocates all the edge records specified in {rd},
    and assigning to them sequential edge id numbers starting from the 
    {edge_id_min} parameter found in the file's header (which is returned in {*eid0_P}). The stores 
    into {A.e[0..A.ne-1]} the base arcs of those edges.  It also stores into {R.e[0..R.ne-1]}  
    the root arcs specified in the file.  The vectors {*A} and {*R} are expanded and trimmed
    as needed. */

#endif
