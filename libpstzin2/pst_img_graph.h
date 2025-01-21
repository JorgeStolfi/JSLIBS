/* Last edited on 2025-01-13 08:02:10 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_H
#define pst_img_graph_H

/* Improved version of {pst_graph.h} that uses {haf.h} for the topology info. */ 

#include <haf.h>
#include <r2.h>
#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_path.h>

#define DEFAULT_WMAG 1.0
  
typedef enum
  { pst_img_graph_mark_NONE,
    pst_img_graph_mark_TO_DELETE,
    pst_img_graph_mark_TO_KEEP,
    pst_img_graph_mark_DELETED
  } pst_img_graph_mark_t;
  /* A vertex or edege state mark. */
  
#define DELETED pst_img_graph_mark_DELETED
  /* Temporary short name. */

typedef struct pst_img_graph_vertex_data_t
  { haf_arc_t aout;              /* One of the directed edges out of the vertex. */
    int32_t x;                   /* Orignal col index in height map, or {-1} if none. */
    int32_t y;                   /* Orignal row index in height map, or {-1} if none. */
    r2_t coords;                 /* Plot coordinates. */
    pst_img_graph_mark_t vmark;  /* Used during enumeration of connected components. */
  } pst_img_graph_vertex_data_t;
  /* The {coords} are used only for plotting the graph.
  
    The fields {x,y}, when not {-1}, are indices of an element of a height map.
  
    The {vmark} is used internally when shrinking or enumerating connected parts
    of the the graph.  It also defines the vertex color when plotting. */

typedef struct pst_img_graph_edge_data_t
  { uint32_t org[2];            /* Index to the origin vertex in each direction. */
    double delta;               /* Height difference in the base direction. */
    double weight;              /* Edge weight (for both directions). */
    pst_img_graph_mark_t emark; /* Used for connected component enumeration. */
    char *label;                /* Label for printout. */
    pst_path_t path;            /* Path for plotting the edge in the base direction.  */
  } pst_img_graph_edge_data_t;
  /* The {path} is used for plotting the edges.
  
    The {emark} is normally {NONE} for valid edges, and {DELETED} for
    edges that were deleted from the graph. Other values are used
    internally when scanning the graph, and define the vertex color when
    plotting. */

typedef struct pst_img_graph_t
  { uint32_t NV;       /* Actual number of of vertices (including NULLs). */
    uint32_t NE;       /* Actual number of (undirected) edges (including NULLs). */
    uint32_t NV_max;   /* Allocated space for vertices */
    uint32_t NE_max;   /* Allocated space for edges */
    pst_img_graph_vertex_data_t *vdata; /* Vertex data, indexed {0..NV-1}*/
    pst_img_graph_edge_data_t *edata;  /* Edge data, indexed {0..NE-1} */
    haf_edge_t *hedge; /* Edges of the half-edge mesh, indexed {0..NE-1}*/
  } pst_img_graph_t;
  /* The vertices are {vdata[0..NV-1]}. */

pst_img_graph_t *pst_img_graph_new(uint32_t NV_max, uint32_t NE_max);
  /* Creates a new graph record {g}. Allocates the lists {g.vdata},
    {g.edata}, and {g.hedge} with sizes {NV_max}, {NE_max}, and {NE_max}
    but initializes {g.NV} and {g.NE} to zero. */

void pst_img_graph_free(pst_img_graph_t *g);
  /* Reclaims all storage used by the graph {g}, including the {haf_edge_t}
    records, the vertex and edge data records, the labels, and the
    plot paths. */

uint32_t pst_img_graph_add_vertex
  ( pst_img_graph_t *g,
    int32_t x,
    int32_t y,
    haf_arc_t aout,
    r2_t coords 
  );
  /* Adds to {g} a vertex with height map indices {x,y} and plot coordinates {coords}.
    The arc {aout} must be one arc from the mesh that leaves the vertex. Returns the
    index of the vertex in the graph. */

haf_arc_t pst_img_graph_edge_add
  ( pst_img_graph_t *g,
    uint32_t org,
    uint32_t dst,
    double d,
    double w,
    char* label,
    pst_path_t path
  );

int32_t pst_img_graph_get_edge_id(haf_arc_t e);
  /* Returns the ID of the undirected edge of the mesh 
    underlying the arc {e}. If {a} is {NULL}, returns {-1}. */

int32_t pst_img_graph_get_arc_id(haf_arc_t a);
  /* Returns the ID of the arc {a}, which is twice the
    ID of the undirected edge {e} plus the longitudinal 
    orientation bit of {a}. If {a} is {NULL}, returns {-1}. */

uint32_t pst_img_graph_get_arc_origin(pst_img_graph_t *g, haf_arc_t a);
  /* Returns the vertex that is the origin of the arc {a}
    (which must not be {NULL}. */

void pst_img_graph_set_arc_origin(pst_img_graph_t *g, haf_arc_t a, uint32_t org);
  /* Sets the vertex with index {org} as the origin of the arc {a}
    (which must not be {NULL}. */

double  pst_img_graph_get_edge_weight(pst_img_graph_t *g, haf_arc_t a);
  /* Gets the weight of the (undirected) edge underlying arc {a}. */

void pst_img_graph_set_edge_weight(pst_img_graph_t *g, haf_arc_t a, double w);
  /* Sets the weight of the (undirected) edge underlying arc {a} to {w}. */

double pst_img_graph_get_arc_delta(pst_img_graph_t *g, haf_arc_t a);
  /* Gets the supposed height difference (/delta/) along the (directed) arc {a}.
    The result for {haf_sym(a)} will be the negative of this. */

void pst_img_graph_set_arc_delta(pst_img_graph_t *g, haf_arc_t a, double delta);
  /* Sets {delta} as the supposed height difference along the 
    (directed) arc {a}. Also sets {-delta} as the difference along {haf_sym(a)}. */

pst_path_t pst_img_graph_get_edge_path(pst_img_graph_t *g, haf_arc_t a);
  /* Gets the plotting path associated with the directed arc {a}. */

void pst_img_graph_set_edge_path(pst_img_graph_t *g,haf_arc_t a, pst_path_t p);
  /* Sets {p} as the plotting path associated with the directed arc {a}.
    Also creates the reversed path and associates it with the 
    arc {haf_sym(a)}. */

pst_path_t pst_img_graph_compute_star_wedge_path
  ( pst_img_graph_t *g,
    uint32_t vi0,
    uint32_t vi1,
    uint32_t vi2
  );
  /* Computes a {pst_path_t} that is suitable for any of the edges 
    that will result from a star-wedge (star-delta) transformation
    involving the three corner vertices {vi0,vi1,vi2}. It causes the edge
    to curve towards the barycenter of the three vertices. */

uint32_t pst_img_graph_vertex_out_degree(pst_img_graph_t *g, uint32_t vi);
  /* Returns the number of arcs out of {g.vdata[vi]}. */

haf_arc_t pst_img_graph_get_connecting_arc(pst_img_graph_t *g, uint32_t vi0, uint32_t vi1);
  /* Returns the arc that goes from vertex number {vi0} to vertex number {vi1},
    of {NULL} if there is no such arc. */

haf_arc_t pst_img_graph_find_rightmost_arc(pst_img_graph_t *g, haf_arc_t a);
  /* Finds the arc out of {a}'s origin that goes to the vertex with
    the highest {X} coordinate. */

void pst_img_graph_check_consistency(pst_img_graph_t *g);
  /* Runs some validity checks on the graph {g}.  Bombs out with messages
    if they fail. */
  
bool_t pst_img_graph_equal(pst_img_graph_t *g, pst_img_graph_t *h);
  /* Compares the two graphs {g} and {h}. If they are isomprphic and 
    isometric, returns {TRUE}.  If there are substantial
    differences, writes them to {stderr} and returns {FALSE}. */
  
uint32_t pst_img_graph_find_nearest_vertex(pst_img_graph_t *g, r2_t *p);
  /* Returns the index of the valid (not deleted) vertex of {g} that is closest to {p}.
    Fails if {g} has no valid vertices. */

double pst_img_graph_compute_left_face_curl(pst_img_graph_t *g, haf_arc_t a);
  /* Estimates the net integral of the curl of the gradient field inside 
    the left face of arc {a}, as the net sum of the deltas of the arcs
    in that face. */

r2_t pst_img_graph_left_face_baricenter(pst_img_graph_t *g, haf_arc_t a);
  /* Computes the barycenter of the left face of arc {a}, as the simple average of 
    the coordinates of its vertices.  Currently ignores the edge paths. */

void pst_img_graph_put_curl_into_image(pst_img_graph_t *g, float_image_t* U);
  /* Stores into the image {U} the curl of the faces of {g}. */

#undef DELETED

#endif


