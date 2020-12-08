#ifndef gem_H
#define gem_H
/* Last edited on 2015-12-01 15:10:51 by stolfilocal */

#define gem_H_copyright "Copyright © 2014 State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <vec.h>

/*
  THE GEM DATA STRUCTURE FOR {d}-DIMENSIONAL COLORED TRIANGULATIONS

    For further details, see the technical report IC-06-016, "Gems: A general
    data structure for {d}-dimensional triangulations" by Arnaldo
    Montagner and Jorge Stolfi, Institute of Computing, UNICAMP,
    September 2006.
    
  SIMPLICIAL MAPS
    
    For any natural number {d}, a /{d}-dimensional simplicial map/ or
    /{d}-triangulation/ is a partition {T} of a compact topological
    space {S} into finite collections of /{k}-elements/, subsets of {S}
    which are homeomorphic to {k}-dimensional open balls, for each {k}
    in {0..d}; with certain incidence constraints between them. 
    
    The 0-elements and 1-elements are called /vertices/ and /edges/,
    respectively. A vertex is a set containing a single point of {S}.
    The {d}-elements are called /cells/, and the {d-1}-elements are
    called /walls/.
    
    The map constrains require that the closure {K*} in {S} of any
    {k}-element {K} of a triangulation {T} be the union of elements of
    {T} (the /faces/ of {K}); and, moreover, the partition of {K*}
    consisting of those elements must be /simplicial/,
    that is, homeomorphic to the geometric facial
    partition of the /closed canonical {k}-simplex/. The latter is the
    convex hull of the points of {R^{k+1}} that comprise its canonical
    basis:
    
    | {(1,0,...,0,0)  (0,1,...,0,0), ..., (0,0,...,1,0), (0,0,...,0,1)
    
    In particular, the boundary of an edge in {S} must be two distinct
    vertices, its /endpoints/; and, for any {k > 0}, the boundary in {S}
    of any {k}-dimensional element {K} contains {k+1} distinct elements of dimension
    {k-1} (the /facets/ of {K}) and {k+1} distinct vertices.
    
    The map constraints also require that every {k}-element with {k < d}
    must be in the boundary of some cell. (It follows that, if the
    triangulation has any elements, it must have at elast one cell.)
    And, finally, that if a {k}-element {E} is in the boundary of two
    distinct cells {A,B}, any partition of the cells of {T} that
    separates {A} and {B} must separate two cells {A1,B1} that have a
    common wall incident to {E}.
    
    It follows from the map constraints that every point {p} that lies on an
    element of dimension {d-2}, {d-1}, or {d} has a neighborhod that is 
    homeomorphic to either the unit {d}-dimensional open ball, or to the
    intersection of that ball with a closed half-space that touches the origin;
    with {p} mapped to the origin, in either case.  Thus the union of all
    elements with those dimensions is a {d}-dimensional manifold with border
    contained in {S}.
    
    A wall is said to be a /border/ wall (or or /free/, /unshared/ or
    /unattached/) if it is incident to only one cell, otherwise it is
    said to be /internal/ (or /shared/, or /attached/). In general, an
    element that is not a cell is said to be a /border/ element if it
    lies on the boundary of a border wall, and /internal/ otherwise. The
    cells are always internal elements.
    
    In particular, a /0-triangulation/ has only a finite set of vertices.
    
    As a special case one also defines a /{-1}-triangulation/ as being
    an empty set of elements, which comprise a partition of the empty
    topological space.
    
    The /standard extension/ creates a {d+1}-triangulation {T'} (unique
    up to homeomorphism) from a {d}-triangulation {T}, by adding one
    vertex {V} to each connected component of {T}, and adding for each
    {k}-element {K} of that component a {k+1}-element {K'} that is the
    open cone of {K} with {V}.
    
    In the standard extension {T'}, the connected components of {T'} are
    in one-to-one correspondence with those of {T}. The incidence
    relations of the new cone elements are the same as those of the
    corresponding original elements, and every one of them is incident
    to the new vertex {V} of its connected component. The border
    elements of {T'} will be all elements of {T}, plus the cones of
    border elements of {T}, plus the new apex {V} of any component that
    has border elements. In particular, there will be a 1-1
    correspondence between the cells of {T} and those of {T'}, and each
    cell of {T} will be a border wall of the corresponding cell of {T'}.
    
    The standard extension can be applied repeatedy to turn any
    {d}-triangulation {T} into a {d+k}-triangulation {T'}, for any
    natural number {k}. In {T'}, for each component of {T} there will be
    a new distinct {k}-face {U} whose closure {U*} is not incident to
    any element of {T}. Every element of {T'} will be an element {K} of
    {T}, a face {F} of {U}, or a generalized cone between some such {K}
    and some such {F} (whose closure is simplicial).
    
    The /empty triangulation/ consists of the only partition of the
    empty topological space, namely the empty set (a partition with zero
    parts). The empty triangulation is a {d}-dimensional triangulation
    for any natural number {d}; it is the standard {d}-extension of
    itself. It is convenient also to define the (only) {d}-triangulation
    with {d == -1} as being the empty triangulation.
    
  COLORED SIMPLICIAL MAPS
    
    A {d}-triangulation is /colored/ if each vertex has an associated
    integer in {0..d} (its /vertex-color/), such that the endpoints of
    every edge have distinct colors; or, equivalently, such that the
    vertices on the boundary of any cell have {d+1} distinct colors.
    
    In such a map, each facet (wall) of each cell also has a /wall-color/ which
    is, by definition, the color of the opposite vertex. The {d+1} walls
    in the boundary of any cell have different wall-colors. The
    /{i}-wall/ of a cell is the wall with wall-color {i}.  The 
    wall-color of a wall is the only vertex-color that does not
    appear on the vertices of that wall.
    
    Note that in a colored 1-triangulation the cells are topological line
    segments, whose endpoints have vertex-colors 0 and 1; and the end
    with wall-color 0 is the vertex with vertex-color 1, and
    vice-versa.
    
    A colored 0-triangulation (where each cell has a single corner) 
    and the empty {-1}-triangulation (that has no cells and no vertices)
    are colored, vacuously.
    
    The standard extension allows one to see a colored {d}-triangulation
    {T}, for any {d>=-1}, as a colored {d+1}-triangulation {T'}. Namely,
    vertex-color {d+1} is assigned to the vertices of {T'} that are not
    vertices of {T}. (one such vertex in each connected component of
    {T'}), and wall-color {d+1} assigned to the cells of {T} (that
    become border walls of the cells of {T'}). */

typedef struct gem_node_t* gem_ref_t;
  /* A value of type {gem_ref_t} is ether {NULL}, or a pointer to a node
    of a /{d}-dimensional GEM data structure/, or /{d}-gem/, for some 
    natural number {d}.
    
    Each node {p} of a {d}-gem has {d+1} fields that are non-null {gem_ref_t}
    pointers, its /adjacency links/. The only constraint on those links
    is the /gem invariant/: if field {i} of a node {p} points to a 
    node {q}, then field {i} of {q} points back to {p}.  
    
    Any link in a node {p} may be a /loop/ that points back to {p} (and
    thus automatically satisfies the gem invariant).
    
    A /trivial {d}-gem/ has only one node, with all its link pointing to
    itself. 
    
    If {d} is positive, a {d}-gem {G} represents some colored
    {d}-triangulation {T}.  Each node {p} of {G} corresponds to a cell
    {K} of {T}. If a link with index {i} of {p} points to a distinct
    node {q}, that node represents a cell {L} of {T} (the only one) that
    shares with {K} the wall with wall-color {i}. If link {i} of {p} is
    a loop, the wall with wall-color {i} of {K} is unshared (in the
    border).
    
    The case {d=0} is somewhat special. Each node of a 0-gem has only
    one adjacency link, with index 0. If that link is always a loop (that is, every
    node of {G} is a trivial 0-gem), the gem {G} can be
    interpreted as a 0-triangulation, that is, a finite set of
    isolated vertices.  If the 0-link of some node is not a loop, that
    node and its neighbor constitute a two-node connected 0-gem,
    which has no corresponding 0-triangulation.  Such a structure will be
    called a /bivertex 0-gem/.
    
    With these definitions, the trivial {d}-gem represents the trivial
    {d}-triangulation, for every {d >= 0}.
    
    The /empty gem/ is the empty set of nodes; it is vacuously a {d}-gem
    for all {d>=0}. It is convenient to define the (only) {d}-gem with
    {d=-1} as being the empty gem, and let it represent the empty
    triangulation (which is a colored {d}-triangulation for all
    {d>=-1}).
    
    A {d}-gem {G} that represents some {d}-triangulation {T} can also be
    viewed as a {d'}-gem of its {d'}-dimensional extension {T'}, for any
    {d'>=d}. In that case the adjacency links with indices {d+1..d'} in
    all nodes are assumed to be all loops. Conversely, if the links with indices
    {d'+1..d} of a {d}-gem are all loops, that gem can be viewed as a {d'}-gem
    of a {d'}-triangulation {T'}, whose {d}-extension is {T}. This is
    true for all {d>=-1}, unless {d>0}, {d' == 0}, and
    
    Note that any {d}-extension of the trivial 0-gem is a trivial {d}-gem;
    and any {d}-extension of the empty gem is the empty gem.
    */
    
#define gem_DIM_MAX 120 
  /* The maximum dimension of a gem, just for consistency checks. */

gem_ref_t gem_node_new(int d);
  /* Creates a trivial {d}-gem consisting of a new single node {p}
    with all the links pointing to {p} itself.  It represents
    the trivial colored triangulation of any non-negative dimension.  
    
    The parameter {d} must be in {0..gem_DIM_MAX}. Only the first {d+1}
    links are actually stored in the node and can be modified (see
    {gem_splice} below). */
   
int gem_node_dim(gem_ref_t p);
  /* Returns the allocated dimension of node {p}, that is, the argument
    {d} that was given to {gem_node_new} when {p} was created. If {p} is
    {null}, returns {-1}. */
 
void gem_node_free(gem_ref_t p);
  /* Detaches the {d}-dimensional gem node {p} from the rest of the structure,
    and reclaims its space.  Any walls of other nodes that were shared with {p} will
    become part of the free border.  If {p} is null, does nothing. */
    
void gem_component_free(gem_ref_t p);
  /* Reclaims the space allocated for all nodes in the entire connected
    component of the gem that owns node {p}. All the nodes will be
    detached from all neighbors before being freed. If {p} is null, does nothing.*/
    
gem_ref_t gem_step(gem_ref_t p, int i);
  /* Walks through the adjacency link {i} of node {p}, returns the
    destination node (which may be {p} itself). 
    
    The index {i} must be non-negative, and {p} must be non-null. If {i}
    is greater than the parameter {d} given to the call of
    {gem_node_new} that created {p}, the link is assumed to be a loop
    and the result is {p} itself.
    
    In the triangulation, this operation moves from the cell represented
    by {p} across its {i}-wall into the adjacent cell {q}, if that wall
    is interior; or stays put, if that wall is in the border.
    
    In any case, this operation is its own inverse:
    {gem_step(gem_step(p,i),i) = p} for any node {p} and any
    non-negative {i}. */

/* EXTENSIBE VECTORS OF NODE REFS */

vec_typedef(gem_ref_vec_t,gem_ref_vec,gem_ref_t);
  /* Defines the type {gem_ref_vec_t}. A variable {v} of that type is a
    descriptor containing a pointer {v.e} to a vector of {gem_ref_t}
    elements, and the allocation size {v.ne} of that vector. 
    
    Use {v = gem_ref_vec_new(N)} to allocate a vector initially with {N}
    elements (that is, with {v.ne == N}. 
    
    Use {gem_ref_vec_expand(&v,k)} to make sure that {*(v.e)} has the
    element with index {k} (that is, to ensure {v.ne > k}), reallocating
    and copying {*(v.e)} if necessary. 
    
    Use {gem_ref_vec_trim(&v,N)} to truncate or expand {*(v.e)} to
    exactly {N} elements (that is, to make {v.ne == N}) preserving its
    first {min(v.ne,N)} elements.
    
    Note that when {gem_ref_vec_expand(&v,k)} reallocates {*(v.e)} it
    will roughly double its size. Therefore, a sequence of {N} calls,
    starting with {v.ne == 0} and with {k} increasing by 1 each time,
    will take time proportional to {N} (or constant amortized time per call),
    not n {N^2}; but may use 3 times as much memory as the minimum needed
    while resizing, and 2 times as much while .
    
    The procedure {gem_ref_vec_trim(&v,N)} however has cost proportional
    to {N} (unless {v.ne} is {N} already), and therefore should be used
    sparingly. */

/* DATA */

void gem_set_data(gem_ref_t p, int data);
int gem_get_data(gem_ref_t p);
  /* Sets and gets the integer {data} field of node {p}. */

/* SPLICING */

void gem_splice(gem_ref_t a, gem_ref_t b, int i);
  /* Joins or separates two nodes {a} and {b} by their {i}-links.  The
    nodes must be non-null, and the index {i} must be in {0..d} where {d
    = min(gem_node_dim(a),gem_node_dim(b))}.
  
    If {a} and {b} represent cells {A} and {B} of some colored triangulation 
    {T} this operation joins or separates {A} and {B} at their
    respective {i}-walls.
    
    More precisely, if the {i}-walls of cells {A} and {B} are both in the border
    (that is, {gem_step(a,i) == a} and {gem_step(b,i) == b}), the procedure glues
    the two walls together so that they become a single wall shared by
    both cells.  Otherwise, if {a} and {b} are adjacent across the
    {i}-wall (that is, {gem_step(a,i) == b} and {gem_step(b,i) == a}),
    the procedure separates the two cells, splitting that shared wall
    into two border walls. All other situations are illegal.
    
    This operation is always its own inverse; that is, the sequence
    {gem_splice(a,b,i);gem_splice(a,b,i)} has no effect. It is also
    symmetric, that is, {gem_splice(a,b,i)} and {gem_splice(b,a,i)} have the
    same effect. */

/* TRAVERSAL AND DOMAINS

  Let {R[0..nr-1]} is a list of non-negative integers. An /{R}-domain/
  of a gem is a set of cells that are pairwise reachable by walking
  across walls whose wall-colors are in that list.
  
  If the nodes in a gem were allocated with dimension parameter {d},
  then any color in the list that is not in {0..d} will be ignored.
  
  The procedures in this section build a list {vis[...]} of gem nodes
  that can be reached from some given starting node. These procedures
  use an integer /label/ field in each node to identify nodes that have
  been found, by assigning to the label of each node its sequential
  index in the list {vis}. A node {p} such that {vis[get_label(p)] == p}
  is said to be /visited/. */
    
int gem_get_label(gem_ref_t p);
  /* Returns the {label} field of {p}, which is set by {gem_traverse} and 
    other traversal routines. */

void gem_domain_traverse(gem_ref_t root, int nr, int R[], gem_ref_vec_t *visP, int *nvisP);
  /* Given a list {R[0..nr-1]} of non-negative integers and a node
    {root} of some gem {G}, finds all unvisited nodes in the {R}-domain
    of {G} that has the node {root}, that must be non-null.
    
    If {G} represents a triangulation {T}, finds all unvisited cells of
    {T} that can be reached from the cell represented by {root} by
    walking through zero or more walls with wall-colors in {R[0..nr-1]}.
    
    The procedure stores the the new nodes it found into the vector {vis == *visP} and 
    increments {*nvisP} with the count of visited nodes. Every such node will have
    its label set to its index in the {vis} vector.
    
    More precisely, on entry {*nvisP} must be a non-negative count {np}. 
    The procedure assumes that the nodes {vis.e[0..np-1]}, and all nodes
    reachable from them, were previously visited, and in particular
    each of them already has its label set to its index in
    {vis.e}. 
    
    The {nv} nodes visited by the procedure itself are stored in
    {vis.e[np..np+nv-1]}, expanding the vector {vis.e} if needed. In
    particular, {vis.e[np]} will be {root}, unless {root} is already
    among the visited nodes {vis.e[0..np-1]}. The preocedure does NOT
    trim the vector {vis.e} to {np+nv} elements; so {vis.ne} may be
    greater than {*nvisP} on output).
    
    In particular, if {nr} is zero, the procedure will merely store the
    {root} node in {vis.e[np]}, if it is not visited already, 
    and increment {*nvisP} accordingly.

    If {R} is {NULL}, the procedure assumes {R[i] = i} for {i} in {0..nr-1}.
    
   ??? Change definition to add a {label_offset} parameter. ???
   !!! Change clients. !!! */
    
void gem_traverse(gem_ref_t root, int d, gem_ref_vec_t *visP, int *nvisP);
  /*  Equivalent to {gem_domain_traverse(root,d+1,NULL,visP,nvisP)}.
    That is, finds and visits all unvisited nodes that are reachable from 
    the non-null node {root} by adjacency links with indices in {0..d}. */

void gem_domains_enum
  ( gem_ref_t root, 
    int nr, int R[], 
    int ns, int S[],
    gem_ref_vec_t *visP,
    int *nvisP,
    gem_ref_vec_t *repP,
    int *nrepP
  );
  /* Finds all the unvisited nodes in a gem
    structure that are accessible from the node {root} through zero or
    more links with indices in {T[0..nt-1]}, where {nt=nr+ns} and {T} is
    the concatenation of the lists {R[0..nr-1]} and {S[0..ns-1]}. 
    The {root} must not be null.
    
    The procedure stores into the vector {vis = *visP} all the newly
    found nodes, updating the count of such nodes {*nvisP}. It also
    stores into the vector {rep = *repP} one representative node from
    each {R}-domain visited, and updates the count of representatives
    {*nrepP}.
    
    Apart from collecting the representative nodes, the procedure is
    equivalent to {gem_domain_traverse(root,nt,T,visP,nvisP)}. In this
    traversal, the links with indices {R[0..nr-1]} are followed with
    higher priority, so each {R}-domain is completely visited before any
    node of the next {R}-domain is visited. In particular, the first
    {R}-domain visited, if any, is the one that has {root}.
    
    Specifically, on entry, {*nrepP} must be a non-negative count {nq}. 
    The procedure assumes that positions {rep.e[0..nq-1]} are already 
    in use, and stores the {nw} representative nodes that it finds 
    itself in {rep.e[nq..nq+nw-1]}, expanding the vector {rep.e} if needed. */

#endif
