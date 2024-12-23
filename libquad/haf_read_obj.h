#ifndef haf_read_obj_H
#define haf_read_obj_H

/* Reading the half-edge representation of 2D meshes. */ 
/* Last edited on 2024-12-22 10:28:01 by stolfi */

#define half_read_obj_H_copyright \
  "Copyright (C) 2023 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <stdio.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <r3.h>

#include <haf.h>

void haf_read_obj_file
  ( FILE *rd,
    haf_edge_id_t eid0, 
    uint32_t *nf_P,
    r3_vec_t *V_P,
    string_vec_t *VL_P,
    haf_arc_vec_t *A_P,
    int32_t **fleft_P,
    int32_t **vorg_P,
    bool_t verbose
  );
  /* Reads from {rd} a description of the surface mesh of a solid object in the OBJ format,
    and and builds the corresponding half-edge structure in memory.
    
    The procedure returns in {*V_P} a vector with the mesh vertex
    coordinates as read from the file namely {V.e[0..nv-1]}
    where {nv = V.ne}. 
    Vertices from the OBJ file that are not used in any face are 
    omitted, so the index of a vertex in {V.e} is not related
    to its index in the file.  
    
    The procedure also stores into {*VL_P} a vector {VL.e[0..nv-1]} of
    vertex labels. Namely {VL.e[k]} will be the label of vertex {V.e[k]},
    read from the '#'-comment at the end of the corresponding "v" line;
    or {NULL} if there was no such comment.
    
    The procedure then infers the edges of the mesh from the face lines
    in the OBJ file. It creates an {haf_edge_rec_t} record {e} for every
    unordered edge that appeara as the side of each face specified in
    the file, and returns in {*A_P} a vector with the base arcs of those
    edges, namely {A.e[0..A.ne-1]}. The edge IDs are assigned
    sequentially starting from {eid0}, and the index of the base
    arc in {A} is the edge's ID minus {eid0}. 
    
    The procedure also stores in {*vorg_P} and {*fleft_P} two newly
    allocated arrays of integers {vorg,fleft}, both with size {na =
    2*A.ne}, that specify the index into {V} of the origin vertex {vorg[ka]}
    and the index of the left face {fleft[ka]} of the arc with ID {ka + 2*eid0}.
    The index of a face is its sequential index in  the OBJ file, 
    counting from 0.
    
    The {OBJ} file may have free borders, consisting of edges that
    are incident to only one face.  Then only one of the two arcs of
    such edge will have a left face.  The {.next} link of the 
    symmetric arc will be {NULL}, and the {fleft} entry will be {-1}.
    
    The procedure aborts with a message if the same oriented arc is used twice
    in the OBJ file.  In particular, this happens if an edge is
    used three or more times, or of two adjacent faces have inconsistent
    orientations, or if a face uses the same edge twice in the same
    orientation.
    
    The vector descriptors {*V_P}, {*VL_P}, and {*A_P} will be overwritten,
    so those vectors should be empty ({.ne = 0}) or unallocated on entry. */

#endif
