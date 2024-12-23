/* Last edited on 2024-12-23 05:50:21 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_img_graph_H
#define pst_img_graph_H

/* Improved version of {pst_graph.h} that uses {oct.h} for the topology info. */ 

#include <oct.h>
#include <r2.h>
#include <float_image.h>
#include <pst_imgsys.h>

typedef struct pst_path_t
  { uint32_t n; /* number of internal vertices*/
    r2_t* v; /*coordinates of internal vertices*/
    bool_t reverse; /* TRUE if the vertices are stored in reversed manner*/
  } pst_path_t;

#define DEFAULT_WMAG 1.0

struct pst_vertex_data_t
  { uint32_t id; /*Unique id*/
    uint32_t mark; /* Used to mark state of the vertex (0) - nothing, (1) - to be removed  (2) - to be preserved */
    r2_t coords; /* Plot coordinates */
    oct_arc_t edge; /*One of the directed edges out of the vertex */
  };

typedef struct pst_vertex_data_t pst_vertex_data_t;

struct pst_edge_data_t
  { uint32_t id;  /*Id - unique*/
    uint32_t org[2];  /*Index to the vextex in the list*/
    uint32_t dst[2]; /*Index to the vextex in the list*/
    double delta[2]; /*Derivative*/
    double weight; /* weight*/
    uint32_t mark; /*marking number*/
    char* label;
    pst_path_t path[2];
  };

typedef struct pst_edge_data_t pst_edge_data_t;

typedef struct pst_edge_t
  { oct_arc_t edge; /*oct edge structure*/
    pst_edge_data_t* data;
  } pst_edge_t;

struct pst_img_graph_t
  { uint32_t n,m; /*actual number of of vertices and edges (including NULLs) */
    uint32_t max_n,max_m; /* Allocated space for vertices and edges */
    uint32_t n_valid,m_valid;
    pst_vertex_data_t* vertex; /*Vertices (or NULLs) (n)*/
    pst_edge_t* edge; /* One directed edge on each undirected edge (m)*/
  };

typedef struct pst_img_graph_t pst_img_graph_t;

void pst_edge_data_free(pst_edge_data_t* ed);

pst_edge_data_t*  pst_edge_data_create(void);

uint32_t pst_img_graph_get_vertex_index_from_image_indices(int32_t ix, int32_t iy, uint32_t NX, uint32_t NY);

void pst_img_graph_get_vertex_image_indices(r2_t *p ,uint32_t NX, uint32_t NY, int32_t *ix, int32_t *iy);

uint32_t  pst_img_graph_vertex_add(pst_img_graph_t* g, uint32_t id, oct_arc_t edge, r2_t coords );

uint32_t pst_img_graph_find_nearest_vertex(pst_img_graph_t* g, r2_t p);

uint32_t pst_img_graph_get_dir_edge_num(oct_arc_t e);

uint32_t pst_img_graph_get_edge_num(oct_arc_t e);

uint32_t pst_img_graph_get_edge_origin(pst_img_graph_t* g, oct_arc_t e);

void pst_img_graph_set_edge_origin(pst_img_graph_t* g, oct_arc_t e, uint32_t org);

double  pst_img_graph_get_edge_weight(pst_img_graph_t* g, oct_arc_t e);

void pst_img_graph_set_edge_weight(pst_img_graph_t* g, oct_arc_t e, double w);

double pst_img_graph_get_edge_delta(pst_img_graph_t* g, oct_arc_t e);

void pst_img_graph_set_edge_delta(pst_img_graph_t* g, oct_arc_t e, double delta);

pst_path_t pst_img_graph_get_edge_path(pst_img_graph_t* g, oct_arc_t e);

void pst_img_graph_set_edge_path(pst_img_graph_t* g,oct_arc_t e,pst_path_t p);

pst_path_t pst_path_create_empty(void);

pst_path_t pst_path_create_single(r2_t coords);

pst_path_t pst_path_reverse(pst_path_t p);

r2_t pst_path_get_vertex(pst_path_t p, uint32_t i);

pst_path_t pst_path_concatenate(pst_path_t p0, r2_t coords, pst_path_t p1);

pst_path_t pst_img_graph_compute_star_wedge_path(pst_img_graph_t* g, uint32_t vi, uint32_t ve, uint32_t vf);

oct_arc_t pst_img_graph_edge_insert
  ( pst_img_graph_t* g,
    uint32_t org,uint32_t dst,
    double d,
    double w,
    char* label,
    pst_path_t path
  );
 
void pst_img_graph_edge_remove(pst_img_graph_t* g,oct_arc_t e);

void pst_img_graph_connect_vertices(pst_img_graph_t* g);

pst_img_graph_t* pst_img_graph_create(uint32_t max_m,uint32_t max_n);

void pst_img_graph_free(pst_img_graph_t *g);

void pst_img_graph_print_vertex(FILE* arq, pst_vertex_data_t* v);

void pst_img_graph_print_edge(FILE* arq,pst_img_graph_t* g, oct_arc_t e);

void pst_img_graph_print(FILE* arq,pst_img_graph_t* g);

oct_arc_t pst_img_graph_find_leftmost_edge(pst_img_graph_t* g, oct_arc_t e0);

uint32_t pst_img_graph_vertex_count_neighbours(pst_img_graph_t* g, uint32_t vi);

oct_arc_t pst_img_graph_check_neighbourhood(pst_img_graph_t* g,uint32_t vi0, uint32_t vi1);

void pst_img_graph_vertex_remove
  ( pst_img_graph_t* g,
    uint32_t vi,
    double* w_i,
    double wmag,
    bool_t verbose
  );
  
void pst_img_graph_reorganize(pst_img_graph_t* g);

void pst_img_graph_check_consistency(pst_img_graph_t* g);
  
pst_img_graph_t* pst_img_graph_copy(pst_img_graph_t* g);

void pst_img_graph_write_vertex(FILE* arq, pst_vertex_data_t* v);
void pst_img_graph_write_path(FILE* arq, pst_path_t p);
void pst_img_graph_write_edge(FILE* arq, pst_img_graph_t *g, oct_arc_t e);
void pst_img_graph_write(FILE* arq, pst_img_graph_t* g);

void pst_img_graph_read_vertex(FILE* arq, pst_img_graph_t * g);
pst_path_t pst_img_graph_read_path(FILE* arq);
void pst_img_graph_read_edge(FILE* arq, pst_img_graph_t *g, int32_t list_onext[], int32_t list_long_bit[]);
pst_img_graph_t* pst_img_graph_read(FILE* arq);

double pst_img_graph_compute_left_face_curl(pst_img_graph_t* g, oct_arc_t e0);
r2_t pst_img_graph_left_face_baricenter(pst_img_graph_t* g, oct_arc_t e0);
void pst_img_graph_put_curl_into_image(pst_img_graph_t* g,float_image_t* OZ);

void debug_vertex_remove_edge(pst_img_graph_t* g, oct_arc_t e);

#endif
