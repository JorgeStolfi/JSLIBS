#ifndef obj_file_read_H
#define obj_file_read_H

/* Basic parsing of Wavefront OBJ mesh files. */ 
/* Last edited on 2025-01-09 23:24:20 by stolfi */

#define obj_file_read_H_copyright \
  "Copyright (C) 2024 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <stdio.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <r3.h>

#include <obj_file.h>

obj_file_data_t *obj_file_read(FILE *rd, bool_t verbose);
  /* Reads from {rd} a solid model in the Wavefront OBJ format.
    Returns the essential information as a {obj_file_data_t} structure {D}.
    
    The procedure returns the mesh vertex coordinates and labels as read
    from the file (from the "v" lines) in the vectors {D.V.e[0..NV-1]}
    and {D.VL.e[0..NV-1]}, the corner texpoint (texture mapping)
    coordinates (from the "vt" lines) in {D.T.e[0..NT-1]}, and the
    corner normals (from the "vn" lines) in {D.N.e[0..NN-1]}. Note that
    indices into these arrays are 1 less than the indices of the things
    in the OBJ file, since the latter starts counting at 1. That is,
    {D.V.e[i]} are actually the vertex numbered {i+1} in the file. These
    vectors will be trimmed so that one will have {NV==D.V.ne==D.VL.ne},
    {NT = D.T.ne}, and {NN = D.N.ne}.
    
    The procedure also returns the face data ("f" lines) into the arrays
    {D.FV}, {D.FT}, and {D.FN}.  The number of faces {NF} will be 
    {D.FV.ne = D.FT.ne = D.FN.ne}
    
    Lines "vt" that have only two texpoint coordinates will get an extra 
    zero coordinate appended to them.
    
    The procedure ignores lines that begin with "o" (object name), "g"
    (grouping), "p" (point), "l" (line), "c_interp", "d_interp",
    "usemtl", or "mtllib". It fails if the file has any other line
    types, including those that specify curved elements (like "vp",
    "curv", "surf", bmat", ...) and other extensions (like "call" and
    "csh").
    
    Currently this function understands only the Unix/Linux end-of-line
    code ("\012" = "\n"). It does not support MS Windows end-of-line
    codes ("\015\012") or Apple ones ("\015"). */
    

#endif
