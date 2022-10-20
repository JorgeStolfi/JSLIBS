/* Types and tools for STL files. */
/* Last edited on 2022-10-20 06:03:02 by stolfi */

#ifndef stmesh_STL_H
#define stmesh_STL_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <i3.h>

#include <stmesh.h>

/* INPUT/OUTPUT */

stmesh_t stmesh_STL_read(char *fileName, bool_t binary, float eps, uint32_t nfGuess, bool_t even, bool_t checkSorted);
  /* Reads the STL file with given name (ascii or binary)
    and converts it to a {stmesh_t} structure. 

    See {stmesh_STL_read_INFO} below for an explanation of how
    the topology of the mesh is recovered from the unstructured STL file.

    The {nfGuess} parameter is a hint for the number of (unoriented)
    faces in the mesh. It is used to pre-allocate tables of faces and
    other items. It can be any number, even zero; however, the procedure
    is more efficient if {nfGuess} is equal to the number of faces, or
    slightly higher.
    
    If {even} is true, each vertex coordinate is quantized by rounding
    to the nearest *even* multiple of the fundamental length {eps}. If {even} is
    false, it is rounded to the nearest integer multiple of {eps}, even or odd.

    if {checkSorted} is true, checks whether the triangles in {mesh.f}
    are sorted by their {.minZ} field in non-decreasing order. Aborts
    with message if not. */
    
#define stmesh_STL_read_INFO \
  "The topology is determined by first rounding every coordinate of every vertex to" \
  " an integer multiple of some length unit {EPS} specified by the user.  Triangle" \
  " corners that become coincident by such rounding are assumed to be the same vertex" \
  " of the mesh.  Any triangle that has one or two coincident vertices is flagged" \
  " and discarded.\n" \
  "\n" \
  "  Then, any two sides of two different triangles that have the same" \
  " vertices as endpoints are assumed to be the same edge of the mesh.  Every edge" \
  " then must be incident to at least one face, usually two, possibly three or more.  Every" \
  " vertex too must be incident to at least two distinct edges and at least one triangle.\n" \
  "\n" \
  "  It is a fatal error if two triangles of the mesh have the same three" \
  " vertices after quantization. Reducing the length unit {EPS} may get" \
  " around such problem; otherwise, the input STL file should be edited" \
  " by removing *both* ofending triangles."

/* GENERIC STL READING */

typedef struct stmesh_STL_r3_t { float c[3]; } stmesh_STL_r3_t;
  /* A point of {\RR^3}, with single-precision coordinates, as in and STL file. */

typedef struct stmesh_STL_face_t 
  { stmesh_STL_r3_t v[3];    /* Triangle vertices. */
    stmesh_STL_r3_t normal;  /* Triangle normal. */
  } stmesh_STL_face_t;
  /* A triangular face of the mesh, as represented in an STL file. */  

typedef void stmesh_STL_face_proc_t(int32_t line, stmesh_STL_face_t *face);
  /* Type of a procedure that is called by {stmesh_STL_read} to process
    each triangle {*face} read from the STL file. The {line} argument is
    a line number to be used in error messages. Note: the face record
    {*face} is reused. */

void stmesh_STL_gen_read(char *fileName, bool_t binary, stmesh_STL_face_proc_t *process_face);
  /* Reads the STL file {fileName} and calls {process_face(line,face)} for each
    face {face} read from it.  
    
    If {binary} is false, assumes ASCII STL format; in that case {line}
    is the line number in the file (counting from 1). If {binary} is
    true, assumes binary STL format; in that case, {line} is 1 fr the
    header, 2 for the number of faces, and is incremented by 1 for
    each face read. */

/* DEBUGGING */

void stmesh_STL_print_triangle(FILE *wr, stmesh_STL_face_t *f); 
  /* Prints the STL face {f} to {wr} for debugging purposes. */

/* VERTEX ROUNDING */

i3_t stmesh_STL_round_point(stmesh_STL_r3_t *p, float eps, bool_t even);
  /* Converts a float-valued point {p} of {\RR^3}, as read from an STL
    file, to an integer point {q} of {\RZ^3}. Namely, rounds each coordinate
    {p.c[k]} of {p} to an integer multiple {eps*q.c[k]} of {eps},
    and returns the integer vector {q}.  If {even} is true,
    the coordinates will be rounded to even integers. */

#endif
