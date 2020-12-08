/* Private type definitions for {stmesh.h} */
/* Last edited on 2016-04-18 01:09:27 by stolfilocal */

#ifndef stmesh_rep_H
#define stmesh_rep_H

#define _GNU_SOURCE
#include <stdio.h>

#include <vec.h>
#include <bool.h>
#include <i3.h>
#include <r3.h>

#include <stmesh_STL.h>
#include <stmesh.h>

/* DATA RECORDS */

typedef struct stmesh_rep_t
  {
    float eps;  /* Fundamental unit of length (mm). */ 

    /* Element counts: */
    uint32_t nv;  /* Number of vertices. */
    uint32_t ne;  /* Number of edges. */
    uint32_t nf;  /* Number of faces. */

    /* Allocation sizes: */
    uint32_t nv_max;  /* Max number of vertices. */
    uint32_t ne_max;  /* Max number of edges. */
    uint32_t nf_max;  /* Max number of faces. */

    /* Element storage: */
    struct stmesh_vert_rep_t *v;  /* The vertices, indexed {0..nv-1}. */
    struct stmesh_edge_rep_t *e;  /* The edges, indexed {0..ne-1}. */
    struct stmesh_face_rep_t *f;  /* The faces, indexed {0..nf-1}. */

    i3_t minQ, maxQ;   /* Quantized bounding box corners. */
  } stmesh_rep_t;
  /* The mesh header record.
  
  The fields {.nv_max,.ne_max,.nf_max} are the total number of 
  entries allocated in {.v,.e,.f}.  The actual counts of 
  vertices, edges, and faces are {.nv,.ne,.nf}, respectively. */
   
typedef struct stmesh_vert_rep_t
  { i3_t pos;    /* Coordinates, quantized. */
  } stmesh_vert_rep_t;
  /* A vertex of the mesh is defined by a quantized point of {\RR^3}, with
    coordinates {eps*pos[0..2]}. */

typedef struct stmesh_edge_rep_t
  { stmesh_vert_t endv[2];  /* Endpoints in vertex table, in {0..mesh.nv-1}. */
    uint32_t degree;        /* Number of faces incident to this edge. */
  } stmesh_edge_rep_t;
  /* An UNORIENTED edge of the mesh is defined by two distinct vertices (endpoints),
    {.endv[0..1]}.  In an UNORIENTED edge, 
    the vertices are sorted so that their indices in {mesh.v} are in increasing order.  
    
    A pointer {e} to an ORIENTED edge (a {stmesh_edge_t}) 
    is {ue|d} where {ue} is a pointer to the the corresponding {stmesh_edge_rep_t}
    record, and {d} is the /direction bit/ of {e}. In the /natural/ 
    orientation (direction bit {d = 0}), the edge is traversed from {ue.endv[0]}
    (origin) to {ue.endv[1]} (destination). In the /reversed/ orientation,
    ({d = 1}) it is traversed the other way: the origin of {e} is
    {ue.endv[1]}, and the destination of {e} is {ue.endv[0]}. */

typedef struct stmesh_face_rep_t
  { stmesh_edge_t side[3];   /* Oriented edges along the perimeter. */
    int32_t minZ;            /* Minimum {Z} coordinate of any vertex, as multiple of {eps}. */
    int32_t maxZ;            /* Maximum {Z} coordinate of any vertex, as multiple of {eps} . */
  } stmesh_face_rep_t;
  /* An UNORIENTED face of the mesh is a triangle defined by its three distinct edges 
    and three distinct vertices.  
    
    The fields {.side[0..2]} are the ORIENTED edges along the 
    perimeter of the triangle, in consistent order.  Namely, the destination
    of the oriented edge {.side[k]} is the origin of the oriented edge
    {.side[(k+1)%3]}, for {k} in {0..2}.
    
    In an UNORIENTED face, the cyclic order of the sides and the
    starting side are chosen so that the ORIENTED index of {side[0]} is
    minimized. Since all six oriented edges on the perimeter are
    distinct, that orientation is unique.
    
    A pointer {f} to an ORIENTED face is {uf | 2*k + d},
    where {uf} is a pointer to the the corresponding {stmesh_face_rep_t}
    record, {k} is a /side selector/ in {0..2}, and {d} is a /direction bit/
    that specifies a sense of traversal of the perimeter. 
    In the /natural direction/ ({d = 0}), the base edge of {f} is {uf.side[k]},
    and the sides are traversed in increasing order, starting from that side -- namely,
    {uf.side[k]}, {uf.side[(k+1)%3]}, and {uf.side[(k+2)%3]}. In the /reversed
    orientation/ ({d = 1}), the sides are traversed in the opposite order
    and orientation, also starting from the base edge:
    {EREV(uf.side[k])}, {EREV(uf.side[(k+2)%3])}, {EREV(uf.side[(k+1)%3])},
    where {EREV(x)} is {stmesh_edge_reverse(x)}. */


#endif
