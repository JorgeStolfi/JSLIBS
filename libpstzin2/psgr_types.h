/* Last edited on 2025-04-27 02:16:03 by stolfi */
/* Basic type declarations for {psgr.h}. */

#ifndef psgr_types_H
#define psgr_types_H

#include <stdint.h>

#include <r2.h>
#include <i2.h>
#include <psgr_path.h>

typedef struct psgr_vertex_t { uint32_t ix; } psgr_vertex_t;
  /* Index of a vertex in a graph. */
  
typedef struct psgr_edge_t { uint32_t ix; } psgr_edge_t;
  /* Index of an (unoriented) edge in a graph. */
  
typedef struct psgr_arc_t { uint32_t ix; } psgr_arc_t;
  /* Index of an arc (oriented edge) in a graph. It is {2*e.ix + db} where
    {e} is the the unoriented edge, and {db} is a /direction
    bit/ bit that specifies the direction of the arc. The /base arc/ of
    the edge is the one with {db=0}. */
   
typedef uint8_t psgr_dir_bit_t;
  /* The direction bit of an arc (oriented edge).  Actually, it can be only 0 or 1. */
   
#define psgr_arc_NONE ((psgr_arc_t){ UINT32_MAX })
#define psgr_edge_NONE ((psgr_edge_t){ UINT32_MAX })
#define psgr_vertex_NONE ((psgr_vertex_t){ UINT32_MAX })
  /* Values for {psgr_vertex_t}, {psgr_arc_t}, or {psgr_edge_t}
    that indicates "no such item". */
 
typedef enum
  { psgr_mark_CLEAR,
    psgr_mark_DELETED,
    psgr_mark_KILL,
    psgr_mark_KEEP
  } psgr_mark_t;
  /* A vertex or edge state mark. */
  
typedef struct psgr_t psgr_t;
  /* A structure that describes a planar graph drawn on the plane.
  
    The graph consists of a finite set of /vertices/ and a finite set of /edges/.
  
    Each vertex {v} has an index {v.ix}, assigned consecutively from 0
    when the vertex is added to the graph. The vertices are located at a
    subset of the corners of a rectangular integer grid with lower
    corner {(0,0)} and upper corner {(NX-1,NY-1)} for some positive
    integers {NX} and {NY}. These may be intrepreted as the elements of
    an image (usually a height map) with {NX} cols and {NY} rows, so the
    pixel indices {(x,y)} of each vertex are the col and row indices of the
    correspondng map element.
    
    The properties of each vertex also include a mark that can be set to
    any {psgr_mark_t} value.
    
    Each edge is a line that connects a pair of vertices. Each edge has
    a mark that can be set to any {psgr_mark_t} value. Each edge {e} has
    an edge index {e.ix}, assigned consecutively from 0 as the edge is
    added to the graph.
    
    An /arc/ is an edge with a specific orientation (sense of
    traversal), indicated by a /direction bit/ (either 0 or 1). The
    /symmetric/ of arc is the same edge with the opposite orientation.
    Each arc has an src index which is twice the index of the underlying
    edge plus the direction bit.
    
    An arc has a definite origin vertex and a destination vertex. Other
    properties of an arc include a positive /weight/ and a finite
    /delta/. The symmetric arc has the same weight but its delta has the
    opposite sign. An arc also has a /path/ which is a polygonal line to
    be used when drawing the arc.
    
    The /onext ring/ of a vertex is the circular list of the arcs that
    leave a vertex, sorted counterclockwise by the initial directions of
    their paths.  This list is empty if the vertex has no outgoing arcs. */

#define psgr_MAX_VERTICES (4096*4096)
#define psgr_MAX_EDGES (3*psgr_MAX_VERTICES + 2)
  /* Should be enough. */

#define psgr_MAX_IMG_SIZE 4096
  /* Should be enough. */

#define psgr_FILE_TYPE "psgr_t"
  /* File type for {psgr_write_file} and {psgr_read_file}. */

#define psgr_FILE_VERSION "2025-04-22"
  /* File version for {psgr_write_file} and {psgr_read_file}. */

#endif


