/* Private type definitions for {stpoly.h} */
/* Last edited on 2016-04-18 18:57:47 by stolfilocal */

#ifndef stpoly_rep_H
#define stpoly_rep_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <i2.h>
#include <r2.h>

#include <stpoly_STP.h>
#include <stpoly.h>

/* DATA RECORDS */

typedef struct stpoly_rep_t
  {
    float eps;  /* Fundamental unit of length (mm). */ 

    /* Element counts: */
    uint32_t nv;  /* Number of vertices. */
    uint32_t ne;  /* Number of edges. */

    /* Allocation sizes: */
    uint32_t nv_max;  /* Max number of vertices. */
    uint32_t ne_max;  /* Max number of edges. */

    /* Element storage: */
    struct stpoly_vert_rep_t *v;  /* The vertices, indexed {0..nv-1}. */
    struct stpoly_edge_rep_t *e;  /* The edges, indexed {0..ne-1}. */

    i2_t minQ, maxQ;   /* Quantized bounding box corners. */
  } stpoly_rep_t;
  /* The mesh header record.
  
  The fields {.nv_max,.ne_max,.nf_max} are the total number of 
  entries allocated in {.v,.e,.f}.  The actual counts of 
  vertices, edges, and faces are {.nv,.ne,.nf}, respectively. */
   
typedef struct stpoly_vert_rep_t
  { i2_t pos;           /* Coordinates, quantized. */
    uint32_t degree;    /* Number of edges incident to thei vertex. */
  } stpoly_vert_rep_t;
  /* A vertex of the mesh is defined by a quantized point of {\RR^2}, with
    coordinates {eps*pos[0..1]}. */

typedef struct stpoly_edge_rep_t
  { stpoly_vert_t endv[2];  /* Endpoints in vertex table, in {0..mesh.nv-1}. */
    int32_t minY;            /* Minimum {Y} coordinate of any vertex, as multiple of {eps}. */
    int32_t maxY;            /* Maximum {Y} coordinate of any vertex, as multiple of {eps} . */
  } stpoly_edge_rep_t;
  /* An UNORIENTED edge of the mesh is defined by two distinct vertices (endpoints),
    {.endv[0..1]}.  In an UNORIENTED edge, 
    the vertices are sorted so that their indices in {mesh.v} are in increasing order.  
    
    A pointer {e} to an ORIENTED edge (a {stpoly_edge_t}) 
    is {ue|d} where {ue} is a pointer to the the corresponding {stpoly_edge_rep_t}
    record, and {d} is the /direction bit/ of {e}. In the /natural/ 
    orientation (direction bit {d = 0}), the edge is traversed from {ue.endv[0]}
    (origin) to {ue.endv[1]} (destination). In the /reversed/ orientation,
    ({d = 1}) it is traversed the other way: the origin of {e} is
    {ue.endv[1]}, and the destination of {e} is {ue.endv[0]}. */

#endif
