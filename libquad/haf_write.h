#ifndef haf_write_H
#define haf_write_H

/* Writing the half-edge representation of 2D meshes. */
/* Last edited on 2023-10-05 19:55:11 by stolfi */

#define half_write_H_copyright \
  "Copyright (C) 2023 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <vec.h>

#include <haf.h>

void haf_write_arc(FILE *wr, haf_arc_t a, int32_t width);
  /* Writes the arc {a} to {wr} in the format "{eid}:{db}" where {eid}
    is {haf_edge_id(a)}, and {db} is the direction bit {haf_dir_bit(a)} (0 or 1).
    The {eid} is padded with spaces at left to have at least {width}
    characters. */

#define haf_write_FILE_TYPE "half-edge"
#define haf_write_FILE_VERSION "2023-10-05"

void haf_write_map
  ( FILE *wr,
    haf_edge_count_t ne,
    haf_arc_t a[],
    haf_edge_id_t eid0,
    haf_arc_count_t nr,
    haf_arc_t root[]
  );
  /* Writes to {wr} a description of the half-edge structure.

    The parameter {a[0..ne-1]} must list exactly one arc on each undirected edge.
    The edge id of arc {a[ke]} must be {eid0 + ke}. 

    The arcs {root[0..nr-1]} are selected arcs that are supposed to be 
    handles to the structure.  Typically they comprise eactly one arc
    on each connected component of the structure.
    
    The file format is a header with four lines

      "begin {haf_write_FILE_TYPE} (format of {haf_write_FILE_VERSION})"
      "edges = {ne}"
      "edge_id_min = {eid0}"

    followed by {ne} lines with "{ke} {a[ke].lnext} {a[ke].sym.lnext}",
    for {ke} in {0..ne-1]}, followed by a line
    
       "roots = {nr}"

    then {nr} lines with "{kr} {root[kr]}" for {kr} in {0..nr-1},
    then finally the file trailer line 
   
      "end {haf_write_FILE_TYPE}"
    
  */

#endif
