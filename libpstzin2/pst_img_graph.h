/* Last edited on 2024-12-27 09:28:54 by stolfi */
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

typedef struct pst_img_graph_vertex_data_t
  { int32_t id;      /*Unique id*/
    int32_t mark;    /* Used during enumeration of connected components. */
    r2_t coords;      /* Plot coordinates */
    haf_arc_t edge;   /* One of the directed edges out of the vertex */
  } pst_img_graph_vertex_data_t;
  /* The {coords} are used only for plotting the graph.
    The {mark} has values (0) - nothing, (1) - to be removed  
    (2) - to be preserved. */

typedef struct pst_img_graph_edge_data_t
  { int32_t id;  /*Id - unique*/
    int32_t org[2];  /* Index to the origin vertex in each direction. */
    int32_t dst[2];  /* Index to the destination vertex in each direction. */
    double delta[2];  /* Height difference in each direction. */
    double weight;    /* Edge weight (for both directions). */
    int32_t mark;    /* Used for connected component enumeration. */
    char* label;      /* Label for printout. */
    pst_path_t path[2];
  } pst_img_graph_edge_data_t;
  /* The {path} is used for plotting the edges. */

typedef struct pst_img_graph_edge_t
  { haf_arc_t edge; 
    pst_img_graph_edge_data_t* data;
  } pst_img_graph_edge_t;

typedef struct pst_img_graph_t
  { int32_t NV;       /* Actual number of of vertices (including NULLs). */
    int32_t NE;       /* Actual number of of edges (including NULLs). */
    int32_t NV_max;   /* Allocated space for vertices */
    int32_t NE_max;   /* Allocated space for edges */
    int32_t NV_valid; /* Counts non-{NULL} vertices. */
    int32_t NE_valid; /* Counts non-{NULL} vertices. */
    pst_img_graph_vertex_data_t *vertex; /*Vertices (or NULLs) (NV)*/
    pst_img_graph_edge_t *edge;  /* One directed edge on each undirected edge (NE)*/
  } pst_img_graph_t;
  /* The vertices are {vertex[0..NV-1]}, but some may be {NULL}.
    The number of non-{NULL} vertices is {NV_valid}. 
    Ditto for {edge[0..NE-1]} and {NE_valid}. 
    
    The ID of a vertex {vertex[iv]} is unrelated to its index {iv]}. */

void pst_img_graph_free(pst_img_graph_t *g);

int32_t  pst_img_graph_add_vertex
  ( pst_img_graph_t* g,
    int32_t id,
    haf_arc_t edge,
    r2_t coords 
  );
  /* Adds to {g} a vertex with ID {ìd}, with plot coordinates {coords}.
    The {edge} must be the index of an oriented edge ({haf_arc_t}) 
    in the topology. */

haf_arc_t pst_img_graph_add_edge
  ( pst_img_graph_t* g,
    int32_t org,
    int32_t dst,
    double d,
    double w,
    char* label,
    pst_path_t path
  );

int32_t pst_img_graph_get_edge_num(haf_arc_t e);
  /* Returns the index of the undirected edge underlying the arc {e}. */

int32_t pst_img_graph_get_dir_edge_num(haf_arc_t e);
  /* Returns the index of the arc {e}, which is twice the
    index of the undirected edge plus the longitudinal 
    orientation bit of {e}. */

int32_t pst_img_graph_get_edge_origin(pst_img_graph_t* g, haf_arc_t e);
  /* Returns the vertex that is the origin of the arc {e}. */

void pst_img_graph_set_edge_origin(pst_img_graph_t* g, haf_arc_t e, int32_t org);
  /* Sets the vertex with index {org} as the origin of the arc {e}. */

double  pst_img_graph_get_edge_weight(pst_img_graph_t* g, haf_arc_t e);
  /* Gets the weight of the (undirected) edge underlying arc {e}. */

void pst_img_graph_set_edge_weight(pst_img_graph_t* g, haf_arc_t e, double w);
  /* Sets the weight of the (undirected) edge underlying arc {e} to {w}. */

double pst_img_graph_get_edge_delta(pst_img_graph_t* g, haf_arc_t e);
  /* Gets the supposed height difference (/delta/) along the (directed) arc {e}. */

void pst_img_graph_set_edge_delta(pst_img_graph_t* g, haf_arc_t e, double delta);
  /* Sets {delta} as the supposed height difference along the 
    (directed) arc {e}. */

pst_path_t pst_img_graph_get_edge_path(pst_img_graph_t* g, haf_arc_t e);
  /* Gets the plotting path associated with the directed arc {e}. */

void pst_img_graph_set_edge_path(pst_img_graph_t* g,haf_arc_t e,pst_path_t p);
  /* Sets {o} as the plotting path associated with the directed arc {e}. */

int32_t pst_img_graph_find_nearest_vertex(pst_img_graph_t* g, r2_t p);

pst_path_t pst_img_graph_compute_star_wedge_path
  ( pst_img_graph_t* g,
    int32_t vi,
    int32_t ve,
    int32_t vf
  );
 
void pst_img_graph_edge_remove(pst_img_graph_t* g,haf_arc_t e);

void pst_img_graph_connect_vertices(pst_img_graph_t* g);

pst_img_graph_t* pst_img_graph_create(int32_t NE_max,int32_t NV_max);

void pst_img_graph_print(FILE* wr,pst_img_graph_t* g);

haf_arc_t pst_img_graph_find_leftmost_edge(pst_img_graph_t* g, haf_arc_t e0);

int32_t pst_img_graph_vertex_count_neighbours(pst_img_graph_t* g, int32_t vi);

haf_arc_t pst_img_graph_check_neighbourhood(pst_img_graph_t* g,int32_t vi0, int32_t vi1);

void pst_img_graph_vertex_remove
  ( pst_img_graph_t* g,
    int32_t vi,
    double* w_i,
    double wmag,
    bool_t verbose
  );
  
void pst_img_graph_reorganize(pst_img_graph_t* g);

void pst_img_graph_check_consistency(pst_img_graph_t* g);
  
bool_t pst_img_graph_compare(pst_img_graph_t* g, pst_img_graph_t* h);
  /* Compares the two graphs {g} and {h}. If they are isomprphic and 
    isometric, returns {TRUE}.  If there are substantial
    differences, writes them to {stderr} and returns {FALSE}. */
  
pst_img_graph_t* pst_img_graph_copy(pst_img_graph_t* g);

void pst_img_graph_write(FILE* wr, pst_img_graph_t* g);

pst_img_graph_t* pst_img_graph_read(FILE* rd);

double pst_img_graph_compute_left_face_curl(pst_img_graph_t* g, haf_arc_t e0);
r2_t pst_img_graph_left_face_baricenter(pst_img_graph_t* g, haf_arc_t e0);
void pst_img_graph_put_curl_into_image(pst_img_graph_t* g,float_image_t* OZ);

void debug_vertex_remove_edge(pst_img_graph_t* g, haf_arc_t e);

/* MAPPING VERTICES TO IMAGE PIXELS */

int32_t pst_img_graph_get_vertex_index_from_image_indices(int32_t ix, int32_t iy, int32_t NX, int32_t NY);

void pst_img_graph_get_vertex_image_indices(r2_t *p ,int32_t NX, int32_t NY, int32_t *ix, int32_t *iy);


#endif


