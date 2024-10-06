/* General tools for creating {tosl.h} meshes.  */
/* Last edited on 2024-10-06 17:20:39 by stolfi  */

#ifndef tosl_mesh_H
#define tosl_mesh_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>

typedef struct tosl_mesh_t
  { int32_t NE;             /* Number of unoriented edges in use. */
    tosl_arc_t *Arc;    /* The arcs. */
    int32_t NV;             /* Number of vertices in use. */
    tosl_point_t *Vpos; /* Vertex coordinates. */
    char **Vlab;            /* Vertex labels. */
    int32_t NE_max;         /* Number of unoriented edges allocated. */
    int32_t NV_max;         /* Number of vertices allocated. */
  } tosl_mesh_t;
  /* A mesh with space for {NE} edges and {NV} vertices,]
    with arcs {Arc[0..2*NE-1]}, vertex coordinates {Vpos[0..NV-1]},
    and vertex labels {Vlab[0..NV-1]}.
    
    The vector {Arc} is allocated with {NE_max} entries and 
    the vectors {Vpos,Vlab} are allocated with {NV_max} entries.*/

#define tosl_mesh_MAX_PLANES (8*1024*1024)
  /* Max number of slicing planes allowed. */

#define tosl_mesh_MAX_EDGES (8*1024*1024)
  /* Max number of mesh edges allowed. */

/* EXPLORING */
    
void tosl_mesh_print(FILE *wr, tosl_mesh_t *mesh);
  /* Prints the mesh data to {wr}. */
  
void tosl_mesh_check(tosl_mesh_t *mesh);
  /* Consistency and statistics of the mesh. */

void tosl_mesh_coord_range_get(tosl_mesh_t *mesh, tosl_point_t *vmin_P, tosl_point_t *vmax_P);
  /* Sets each coordinate {*vmin_P} and {*vmax_P} to the min and max of the corresponding
    coordinate among all vertces of the {mesh}. */

void tosl_mesh_coord_range_print(FILE *wr, char *pref, tosl_point_t *vmin, tosl_point_t *vmax, char *suff);
  /* Prints the vertex coordinate ranges {vmin.c[j] .. vmax.c[j]} for {j} in {0..2}
    preceded by {pref} and followed by {suff}. */

void tosl_mesh_arc_print(FILE *wr, char *pref, tosl_arc_id_t ka, char *suff, tosl_mesh_t *mesh);
  /* Prints to {wr} the arc id {ka} and the fields of arc {a =
    mesh.Arc[ka]} The field {mesh.Vpos[a.ivorg]} is printed with
    {tosl_mesh_vert_print}. Everything is preceded by {pref} and
    followed by {suff}. The arc id {ka} and the links {.skip}, {.pred},
    and {.succ} are formatted by {tosl_arc_id_to_string}. */
 
void tosl_mesh_vert_print(FILE *wr, char *pref, tosl_vert_id_t kv, char *suff, tosl_mesh_t *mesh);
  /* Prints to {wr} the vertex id {kv} as "v{kv}" then the label {mesh.Vlab[kv]} (if not {NULL} 
    and coordinates {mesh.Vpos[kv]}. Everything is preceded by {pref} and followed by {suff}. */
 
/* CREATING */

tosl_mesh_t *tosl_mesh_new(int32_t NE_max, int32_t NV_max);
  /* Creates a mesh structure with space for space for {NE_max}
    unoriented edges and {NV_max} vertices. Initially it is empty; that
    is, {mesh.NE = mesh.NV = 0}. */
    
void tosl_mesh_free(tosl_mesh_t *mesh);
  /* Frees all storage associated with the given {mesh}, including the
    internal tables, all vertex labels, and the record {*mesh}
    itself. */

tosl_vert_id_t tosl_mesh_add_vert
  ( tosl_point_t *v,
    char *lab,
    tosl_mesh_t *mesh
  );
  /* Sets {mesh.Vpos[kv] = v} where {kv} is the value of {mesh.NV} on
    entry. If {lab} and {mesh.Vlab} are not null, also sets
    {mesh.Vlab[kv] = lab}. Then increments {mesh.NV}. Returns {kv}.
    Fails if {mesh.NV} is already {mesh.NV_max}. */

tosl_arc_id_t tosl_mesh_add_edge
  ( tosl_vert_id_t kv0,
    tosl_vert_id_t kv1,
    int32_t set_skip,
    tosl_mesh_t *mesh
  );
  /* Sets {mesh.Arc[ia]} to an arc with origin {kv0}, and
    {mesh.Arc[ia+1]} to an arc with origin {kv1}, where {ia} is twice
    the value of {mesh.NE} on entry. The indices {kv0} and {kv1} must be
    distinct and in {0..mesh.NV-1}. 
    
    Also sets the {.skip} fields of those two arcs according to the
    boolean parameter {set_skip}. If {set_skip} is 1 (true), sets the
    links so that the edge becomes an isolated component of the mesh;
    namely, {mesh.Arc[ia].skip = sym(ia)} and {mesh.Arc[sym(ia)].skip =
    ia}. If {set_skip} is 0 (false), sets both {.skip} fields to {-1}.
    
    Also sets the {.succ} and {.pred} so as to make each arc into a
    singleton circular list. Then increments {mesh.NE}. Returns the
    first arc index {ia}. Fails if {mesh.NE} is already
    {mesh.NE_max}. */
    
void tosl_mesh_splice(tosl_arc_id_t ia, tosl_arc_id_t ja, tosl_mesh_t *mesh);
  /* The splice operation for arcs {ia} and {ja}, which must be distinct. */

tosl_arc_id_t tosl_mesh_add_ring
  ( int32_t n,
    char *pref,
    tosl_mesh_t *mesh
  );
  /* Adds to the {mesh} a new connected component that consists of a
    ring of {n} edges, with {n} new vertices and two faces. Returns the
    index of an arc of one of the edges.If {lab} is not {NULL}, the new
    vertices wil have undefined coordinates and labels "{pref}{iv}"
    where {iv} ranges from 0 to {n-1}. Both {mesh.NV} and {mesh.NE} are
    incremented by {n}. Fails if runs out of space in {mesh}. */
    
void tosl_mesh_add_path
  ( int32_t n,
    tosl_arc_id_t ia0,
    tosl_arc_id_t ia1,
    char *pref,
    tosl_mesh_t *mesh
  );
  /* Adds to the {mesh} {n+1} new edges and {n} vertices that connect
    the destination of arc {ia0} to the destination of arc {ia1}, so
    that their {.skip} cycles are broken or joined at those vertices by
    the path.If {lab} is not {NULL}, the new vertices will have labels
    "{pref}{iv}" where {iv} ranges from 0 to {n-1}. The counts {mesh.NV}
    and {mesh.NE} are incremented by {n} and {n+1}, respectively. Fails
    if runs out of space in {mesh}. */
    
void tosl_mesh_link_triang
  ( int32_t parity,
    tosl_arc_id_t ka0,
    tosl_arc_id_t ka1,
    tosl_arc_id_t ka2,
    tosl_mesh_t *mesh
  );
  /* If {parity} is 1, links the arcs {ka0,ka1,ka2}, in that order, as the border of a face,
    through the {.slip} links. If {parity} is 0, links instead {sym(ka2),sym(ka1),sym(ka0)}.
    
    The {.skip} links of the arcs that are to be linked must be {-1}. */

#endif
