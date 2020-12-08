/* Reading and writing structured triangle meshes in STM format. */
/* Last edited on 2016-11-16 15:38:56 by stolfilocal */

#ifndef stmesg_STM_H
#define stmesg_STM_H

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>

#include <stmesh.h>

stmesh_t stmesh_STM_read(FILE *rd);
  /* Reads from the file handle {rd} a structured mesh and
    returns it as an {stmesh_t} structure. See 
    {stmesh_STM_file_fomat_INFO} below for the file format. 
    The file should be opened by the client and will be left open 
    at the end. */

stmesh_t stmesh_STM_read_named(char *fileName, bool_t verbose);
  /* same as {stmesh_STM_read} but reads from  a file named {fileName},
    which is opened an closed by the procedure.  If {verbose}
    is true, writes a message to {stderr}. */
     
void stmesh_STM_write(FILE *wr, stmesh_t mesh);
  /* Writes to the file handle {wr} the given triangle mesh structure
    {mesh}. See {stmesh_STM_file_fomat_INFO} below for the file format. 
    The file should be opened by the client and will be left open 
    at the end. */

void stmesh_STM_write_named(char *fileName, stmesh_t mesh, bool_t verbose);
  /* same as {stmesh_STM_write} but writes to a file named {fileName},
    which is opened an closed by the procedure. If {verbose}
    is true, writes a message to {stderr}. */

#define stmesh_STM_FILE_TYPE "stmesh_t"
#define stmesh_STM_FILE_VERSION "2016-04-30"
  /* Current version of {stmesh} files. */
    
#define stmesh_STM_file_format_INFO \
  "An STMESH (\"structured triangle mesh\") has the general" \
  " structure\n" \
  "\n" \
  "    begin " stmesh_STM_FILE_TYPE " (format of " stmesh_STM_FILE_VERSION ")\n" \
  "    eps = {EPS}\n" \
  "    nv = {NV}\n" \
  "    ne = {NE}\n" \
  "    nf = {NF}\n" \
  "    vertices\n" \
  "    {VERT_DATA}\n" \
  "    edges\n" \
  "    {EDGE_DATA}\n" \
  "    faces\n" \
  "    {FACE_DATA}\n" \
  "    end " stmesh_STM_FILE_TYPE "\n" \
  "\n" \
  "  The line breaks must be as above.  The parameter {EPS} is" \
  " the unit of measure, in millimeters, as a floating-point" \
  " number in '%e' or '%f' format.  The parameters {NV}, {NE}, and" \
  " {NF} are the counts of vertices, edges, and" \
  " faces (triangle), respectively.\n" \
  "\n" \
  "  The {VERT_DATA} consists of one line for each vertex, with" \
  " the vertex index (a sequential number in {0..NV-1}), followed" \
  " by three signed integers {IX}, {IY}, and {IZ} in the" \
  " range {-M..+M} where {M = 2^31-1}, separated by" \
  " blanks.  The coordinates of the vertex, in" \
  " millimeters, are then {(IX*EPS, IY*EPS, IZ*EPS)}.\n" \
  "\n" \
  "  The {EDGE_DATA} consists of one line for each" \
  " unoriented edge, with the edge index (a sequential number" \
  " in {0..NE-1}), followed by two" \
  " integers {UVX[0..1]}, with {UVX[0] < UVX[1]}, separated by" \
  " blanks, that identify the origin and destination of the" \
  " edge, respectively, in its natural orientation.\n" \
  "\n" \
  "  The {FACE_DATA} consists of one line for" \
  " each unoriented face, with the face" \
  " index (a sequential number in {0..NF-1}), followed by" \
  " three undirected edge indices {UXE[0..2]}, with" \
  " {UXE[0] < UXE[1] < UXE[2]}, separated" \
  " by blanks.  The base of the triangle, in its natural" \
  " orientation, will be the edge {UXE[0], in its" \
  " natural orientation.\n" \
  "\n" \
  "  In a well-formed mesh, no two" \
  " triangles should have more than one edge in common, in any" \
  " orientation or order.  No" \
  " two edges should have the same set of endpoints.  No" \
  " two vertices should have the same coordinates."

#endif
