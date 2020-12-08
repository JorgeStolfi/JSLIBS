#ifndef g3map_H
#define g3map_H
/* Last edited on 2016-01-01 18:10:30 by stolfilocal */

/* 
  3-DIMENSIONAL MAPS
  
  A {3}-dimensional map has cells whose closures are homeomorphic to 
  the {3}-dimensional ball, satisfying several conditions.
  
  !!! Describe the conditions !!!
  
  An {3}-dimensional map is decribed internally by
  a gem of dimension {3}, which is the barycentric
  subdivision of the map.  In each simplex of the 
  gem, the corner with color {k} lies in 
  some part of the map with dimension {k},
  for {k} in {0,1,2,3}. The 
  gem may have dimension {d>3}, but in that
  case all nodes share the same corners
  with color {4,5,...d}, and all the walls
  with those colors are unattached.
  
*/

#define _GNU_SOURCE
#include <gem.h>

typedef gem_ref_t g3map_place_t;
  /* A {g3map_place_t} specifies a /place/ on a 3-map. It is actually a
    pointer to one tetrahedron of the barycentric gem, and therefore
    specifies one part of dimension {k}, for each {k} in {0,1,2,3}.
    
    Any place in an 3-map has {4} pointers to other places.
    
    All places on a {k}-part of an 3-map correspond to gem simplices
    with the same corner of color {i}, for all {i >= k}. To reach all
    the places on a {k}-part, one starts from any place on it and
    enumerate all places reached by {i}-pointers with {i < k}. */
    
g3map_place_t g3map_step(g3map_place_t a, int i);
  /* Takes a step on the 3-map from place {a} to the place
    that differs from {a} only in the part of dimension {i},
    which must be in {0..3}.
    
    If {i=0}, moves to the other end of the same edge.
    If {a} is an unattached vertex, returns {a} itself.
    
    If {i=1}, {a} must be a place on some edge (not an
    unattached vertex). Moves to the other edge of the
    same face that has the same origin as {a}.  If that
    endpoint of {a} is  unattached, retuns {a} itself.
    
    If {i=2}, {a} must be a place on some face (not an 
    unattached edge or vertex). Moves to the other face of the same cell
    while remaining near the same edge and near the same
    vertex.  If that edge is a border edge, retrns 
    {a} itself.
    
    If {i=3}, {a} must be a place on some cell (not an
    unattached vertex, edge, or face). Moves to the adjacent cell,
    while remaining near the same face, edge, and vertex. If 
    that face is a border face, returns {a} itself.
    
    The result may be {a} itself only if the map has dimension {i}
    and the part of dimension {i} associated with {a} 
    is a border facet of the map.*/
    
/* MAPS OF DIMENSION 0

  A map of dimension zero is just an isolated vertex. */

g3map_place_t g3map_vert_make(void);
  /* Creates a new 0-dimensional map consisting of a single vertex. The
    repesentation is such that this vertex can be part of any
    3-map. */
    
/* MAPS OF DIMENSION 1

  A map of dimension 1 is a path of one or more edges.
  
  When a path is closed -- that is, every edge has two
  incident edges -- the result is automatically a 
  map of dimension 2, consisting of a single face
  not attached to any cell or any other face. */

g3map_place_t g3map_edge_make(g3map_place_t u, g3map_place_t v);
  /* Creates a new 1-dimensional map consisting of a single edge,
    whose endpoints are the two previously created vertices {u,v},
    which must be distinct and completely unattached.  Returns the
    place nearest to {u}. The representation is such that this
    edge can be part of any 3-map.  */
    
void g3map_edge_splice(g3map_place_t a, g3map_place_t b);
  /* Given two edge places {a,b}, distinct and not attached to any cells,
    will join or separate them by one endpoint.
    
    Specifically, the vertices associated with {a} and {b} may be
    distinct and incident to only one edge each, and then they will be
    identified; or may be the same vertex, shared by two distinct edges,
    and then the edges will be separated at that vertex.
    
    When a cycle of edges is closed, the result is a 2-map -- an
    isolated face with the topology of a closed disk. . If the edges are
    part of the boundary of the same unattached face, that face is
    destroyed and the map becomes unidimensional.
    
    If {a} and {b} are the two ends of the same edge, the 
    result will be a loop edge, enclosing a single face with 
    only one side. */
    
/* MAPS OF DIMENSION 2

  A map of dimension 2 is a surface consisting of one or
  more faces connected by edges, so that there are
  at most two faces incident to each edge. 
  
  When a 2-map becomes closed -- that is, all edges are incident to two
  faces -- the map automatically becomes a 3-map, consisting of a single
  unattached cell whose boundary is that set of faces. The surface
  defined by those faces then had better have the topology of a
  2-sphere, otherwise the cell will have a non-manifold point in its
  interior. */

g3map_place_t g3map_face_make(int n, g3map_place_t e[]);
  /* Creates a new 2-dimensional map consisting of a single 2-face,
    whose sides are the previously created edges {e[0..n-1]},
    which must be all distinct, non-loop, and not attached to any edge, face
    or cell. The origin vertex of edge {e[(i+1)%n]} is identified 
    with the destination vertex of edge {e[i]}. The representation 
    is such that this face can be part of any 3-map.  Returns a place
    on the face near edge end {e[0]}. */
    
void g3map_face_splice(g3map_place_t f, g3map_place_t g);
  /* Given two distinct places {f,g} on two faces, not attached to any cells,
    will join or separate the faces by the edge.
    
    The places {f} and {g} must not be opposite ends of the same edge.
    
    Specifically, the edges associated with {f} and {g} may be distinct
    and incident to only one face each, and then they will be
    identified; or may be the same edge, shared by two distinct faces,
    and then the faces will be separated at that edge.
    
    If the splice creates a closed collection of faces, the result is a
    3-map consisiting of a single unattached cell. Conversely, if the
    two faces were on the boudary of a single unattached cell, that cell
    is destroyed and the result is a 2-map with border. . */
    
/* MAPS OF DIMENSION 3

  A map of dimension 3 is a 3-dimensional quasi-manifold consisting of one or
  more cells connected by shared faces, so that there are
  at most two cells incident to each face. 
  
  When a 3-map becomes closed, nothing magic happens.  There is 
  no 4th dimension. */

void g3map_cell_splice(g3map_place_t a, g3map_place_t b);
  /* Given two distinct places {a,b} on two distinct faces,
    glues or splits two cells by those faces. 
    
    Either the 2-faces of {a} and {b}
    are distinct but isomorphic, and both on the border of the 3-map, 
    or they are the same interior 2-face of the 3-map.  In the first case,
    the two 2-faces are identified, and become a single interior 2-face.
    In the second case, the map is cut apart at the shared 2-face,
    leaving {a} and {b} on two distinct border 2-faces. */

#endif
