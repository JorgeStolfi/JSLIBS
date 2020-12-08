/* Cross-section input/output. */
/* Last edited on 2015-11-16 00:56:03 by stolfilocal */

#ifndef stmesh_section_write_H
#define stmesh_section_write_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <vec.h>

#include <stmesh.h>

typedef struct stmesh_section_t
  {
    float eps;        /* Fundamental unit of length for this section (mm). */
    int32_t pZ;       /* Quantized {Z} corodinate (in multiples of {eps}). */
    uint32_t ns;      /* Total number of line segments in the section. */
    uint32_t nc;      /* Number of separate paths (or loops). */
    r2_t *v;          /* Vertices of segments, grouped by paths. */
    uint32_t *cstart; /* Start and end of each path in {v}. */
  } stmesh_section_t;
  /* Describes a polygonal figure on some horizontal plane, consisting of 
    {ns} straight line segments, partitioned into {nc} polygonal paths,
    open of closed.  The {Z} coordinate of the plane is {pZ*eps}.  The {X} and {Y} 
    coordinates of the vertices are {v[0..nv-1]}, where {nv = ns+nc}.
    It is meant to be used for horizontal slices of the mesh.
    
    The vertices of path {k} are {v[ak..bk-1]} where {ak = cstart[k]} and {bk = cstart[k+1]},
    for {k} in {0..nc-1}.  The path is closed if and only if {v[ak] = v[bk-1]}; in either case,
    the path uses {bk-ak} points of {v} and has {bk-ak-1} segments. 
    Note that {cstart} has {nc+1} elements, with {cstart[0] = 0}
    and {cstart[nc] = nv = ns+nc}.  
    
    The section may have zero paths, but each path has at least one segment
    and two vertices. 
    
    ??? eps should be {double} ??? */

stmesh_section_t* stmesh_section_make(stmesh_t mesh, int32_t pZ, uint32_t nc, uint32_t estart[], stmesh_edge_t e[]);
  /* Builds a cross-section structure from a plane and a set of mesh edges that cross that plane.
    
    The section is assumed to have {nc} open of closed paths. The
    vertices of path {k}, for {k} in {0..nc-1}, are assumed to be the
    intersection of the horizontal plane with quantized {Z}-coordinate
    {pZ} and the edges {e[ak..bk-1]}, where {ak = estart[k]} and 
    {bk = estart[k+1]}.  The path is closed if an only if {e[ak] = e[bk-1]}.
    In any case, the path consists of {nvk = bk - ak} vertices
    and {nvk-1} line segments connecting successive vertices.
    
    Note that {estart} must have {nc+1} elements, with {estart[0]= 0}. 
    The total number of vertices in all paths is {nv = estart[nc]}; the edge list 
    {e} is assumed to have {nv} elements. */

void stmesh_section_free(stmesh_section_t *sec);
  /* Reclaims the space used bu {sec}, including the internal tables 
    and the descriptor record {*sec} itself. */

void stmesh_section_write(FILE *wr, stmesh_section_t *sec);
  /* Writes to file {wr} a cross-section in the format described by the 
    string {stmesh_section_format_INFO} below. */

stmesh_section_t* stmesh_section_read(FILE *rd);
  /* Reads a cross-section from file {rd}, that should be in the format
    described by the string {stmesh_section_format_INFO} below.
    The descriptor and the vertex tables are newly allocated by the
    procedure. */
   
#define stmesh_section_VERSION "2015-11-14"
  /* Current version of cross-section file format. */

#define stmesh_section_format_INFO \
  "The file for each cross-section starts with a preamble\n" \
  "\n" \
  "    begin stmesh_section (format of " stmesh_section_VERSION ")\n" \
  "    eps = {EPS}\n" \
  "    planeZ = {PLANEZ}\n" \
  "    planeQ = {PLANEQ}\n" \
  "    nSegs = {NSEGS}\n" \
  "    nPaths = {NPATHS}\n" \
  "\n" \
  "  and ends with with the file footer line\n" \
  "\n" \
  "    end stmesh_section\n" \
  "\n" \
  "  The parameter {EPS} is a fundamental unit of length (in mm), used to the" \
  " quantize the {Z} coordinate; {PLANEZ} is the {Z}-coordinate" \
  " of the slicing plane (in mm); {PLANEQ} is the integer {PLANEZ/EPS}; and" \
  " {NSEGS} is the total number of line segments in the cross-section.\n" \
  "\n" \
  "  Between the file header and footer there are {NPATHS} (zero or more) blocks of lines, one" \
  " for each separate polygonal path or loop found by the loop closing algorithm.  Each" \
  " block starts with the path header lines\n" \
  "\n" \
  "    path \n" \
  "    nPoints = {NV} \n" \
  "\n" \
  "  and ends with with the path footer line\n" \
  "\n" \
  "    end path\n" \
  "\n" \
  "  In between the path header and footer are the {X} and {Y} coordinates" \
  " of the {NV} path vertices, in millimeters and separated by blanks, one" \
  " vertex per line.  The path can be assumed to be" \
  " closed if the last vertex is equal to the first one.  Note that an" \
  " open or closed path with {NS} sides will have {NV = NS+1} vertices."

#endif
