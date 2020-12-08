/* Last edited on 2013-03-18 01:40:31 by stolfilocal */
/* Created by Rafael F. V. Saracchini */


#include <oct.h>
#include <r2.h>
#include <float_image.h>
#include <pst_imgsys.h>

typedef struct pst_path_t{
  long int n; /* number of internal vertices*/
  r2_t* v; /*coordinates of internal vertices*/
  bool_t reverse; /* TRUE if the vertices are stored in reversed manner*/
} pst_path_t;

#define MARK_VERTEX_NONE 0
#define MARK_VERTEX_REMOVED 1
#define MARK_VERTEX_PRESERVED 2

#define DEFAULT_WMAG 1.0

struct pst_vertex_data_t{
 long int id; /*Unique id*/
 int mark; /* Used to mark state of the vertex (0) - nothing, (1) - to be removed  (2) - to be preserved */
 r2_t coords; /* Plot coordinates */
 oct_arc_t edge; /*One of the directed edges out of the vertex */
};

typedef struct pst_vertex_data_t pst_vertex_data_t;

struct pst_edge_data_t{
  long int id;  /*Id - unique*/
  long int org[2];  /*Index to the vextex in the list*/
  long int dst[2]; /*Index to the vextex in the list*/
  double delta[2]; /*Derivative*/
  double weight; /* weight*/
  int mark; /*marking number*/
  char* label;
  pst_path_t path[2];
  
};

typedef struct pst_edge_data_t pst_edge_data_t;

typedef struct pst_edge_t{
  oct_arc_t edge; /*oct edge structure*/
  pst_edge_data_t* data;
} pst_edge_t;




struct pst_img_graph_t{
  
  long int n,m; /*actual number of of vertices and edges (including NULLs) */
  long int max_n,max_m; /* Allocated space for vertices and edges */
  long int n_valid,m_valid;
  pst_vertex_data_t* vertex; /*Vertices (or NULLs) (n)*/
  pst_edge_t* edge; /* One directed edge on each undirected edge (m)*/
    
};

typedef struct pst_img_graph_t pst_img_graph_t;

void pst_edge_data_free(pst_edge_data_t* ed);
pst_edge_data_t*  pst_edge_data_create(void);

long int pst_img_graph_get_vertex_index_from_image_indices(long int ix, long int iy, long int NX, long int NY);
void pst_img_graph_get_vertex_image_indices(r2_t *p ,long int NX, long int NY, long int *ix, long int *iy);

long int  pst_img_graph_vertex_add(pst_img_graph_t* g, long int id,oct_arc_t edge, r2_t coords );

long int pst_img_graph_find_nearest_vertex(pst_img_graph_t* g, r2_t p);


long int pst_img_graph_get_dir_edge_num(oct_arc_t e);
long int pst_img_graph_get_edge_num(oct_arc_t e);



long int pst_img_graph_get_edge_origin(pst_img_graph_t* g, oct_arc_t e);
void     pst_img_graph_set_edge_origin(pst_img_graph_t* g, oct_arc_t e,long int org);


double  pst_img_graph_get_edge_weight(pst_img_graph_t* g, oct_arc_t e);
void    pst_img_graph_set_edge_weight(pst_img_graph_t* g, oct_arc_t e,double w);


double   pst_img_graph_get_edge_delta(pst_img_graph_t* g, oct_arc_t e);
void     pst_img_graph_set_edge_delta(pst_img_graph_t* g,oct_arc_t e,double delta);


pst_path_t pst_img_graph_get_edge_path(pst_img_graph_t* g, oct_arc_t e);
void pst_img_graph_set_edge_path(pst_img_graph_t* g,oct_arc_t e,pst_path_t p);


pst_path_t pst_path_create_empty(void);
pst_path_t pst_path_create_single(r2_t coords);
pst_path_t pst_path_reverse(pst_path_t p);
r2_t pst_path_get_vertex(pst_path_t p, long int i);
pst_path_t pst_path_concatenate(pst_path_t p0, r2_t coords, pst_path_t p1);
pst_path_t pst_img_graph_compute_star_wedge_path(pst_img_graph_t* g,long int vi, long int ve, long int vf);


oct_arc_t pst_img_graph_edge_insert(pst_img_graph_t* g,
			       long int org,long int dst,
			       double d,
			       double w,
			       char* label,
			       pst_path_t path
 );
 
void pst_img_graph_edge_remove(pst_img_graph_t* g,oct_arc_t e);
void pst_img_graph_connect_vertices(pst_img_graph_t* g);



pst_img_graph_t* pst_img_graph_create(long int max_m,long int max_n);
void pst_img_graph_free(pst_img_graph_t *g);


void pst_img_graph_print_vertex(FILE* arq, pst_vertex_data_t* v);

void pst_img_graph_print_edge(FILE* arq,pst_img_graph_t* g, oct_arc_t e);

void pst_img_graph_print(FILE* arq,pst_img_graph_t* g);

void pst_img_graph_mark_vertex_removal(pst_img_graph_t* g, int degree[]);
void pst_img_graph_mark_vertex_removal_no_bucket(pst_img_graph_t* g);

void pst_img_graph_shrink(pst_img_graph_t* g,double* w_i,double wmag);

oct_arc_t pst_img_graph_find_leftmost_edge(pst_img_graph_t* g, oct_arc_t e0);

long int pst_img_graph_vertex_count_neighbours(pst_img_graph_t* g, long int vi);
oct_arc_t pst_img_graph_check_neighbourhood(pst_img_graph_t* g,long int vi0, long int vi1);
void pst_img_graph_vertex_remove(pst_img_graph_t* g, long int vi,double* w_i,double wmag,bool_t verbose);
void pst_img_graph_reorganize(pst_img_graph_t* g);
void pst_img_graph_remove_paralel_edges(pst_img_graph_t* g);

void pst_img_graph_check_consistency(pst_img_graph_t* g);

pst_imgsys_t* pst_img_graph_build_integration_system(pst_img_graph_t* g,double* iW,long int NX_Z,long int NY_Z,long int** ref_tab);
pst_img_graph_t* pst_img_graph_copy(pst_img_graph_t* g);

void pst_img_graph_solve_system(
    pst_img_graph_t* g,
    pst_imgsys_t* S,
    double *iZ,
    double* iW,
    long int* ref_tab,
    long int maxIter, 
    double convTol, 
    int para, 
    int szero, 
    bool_t verbose 
  );

  void pst_img_graph_copy_solution_from_shrunk(pst_img_graph_t* jg,double* jZ,pst_img_graph_t* ig, double *iZ);
  /* Also marks vertices of {ig} as MARK_VERTEX_PRESERVED or MARK_VERTEX_REMOVED*/
  
  void pst_img_graph_estimate_from_shrunk(pst_img_graph_t* jg,double* jZ,pst_img_graph_t* g, double *iZ);
  
  

void pst_img_graph_integration_recursive( 
  pst_img_graph_t* g,
  double* iZ,
  double* iW,
  double wmag,
  long int maxIter,
  double convTol, 
  int para, 
  int szero, 
  bool_t verbose,
  int level,
  float_image_t* OZ, /*debug only*/
  float_image_t* RZ,
  char* out_prefix,
  bool_t debug
);

void pst_img_graph_integration( 
  pst_img_graph_t* g,
  double* iZ,
  double* iW,
  long int maxIter,
  double convTol, 
  int para, 
  int szero, 
  bool_t verbose,
  int level,
  float_image_t* OZ, /*debug only*/
  float_image_t* RZ,
  char* out_prefix
);

void pst_img_graph_put_solution_into_image(pst_img_graph_t* g, double* iZ, float_image_t* OZ);
void pst_img_graph_put_error_into_image(pst_img_graph_t *g,double* iZ,double* iW, float_image_t *RZ,float_image_t *OZ);

void pst_img_graph_write_vertex(FILE* arq, pst_vertex_data_t* v);
void pst_img_graph_write_path(FILE* arq, pst_path_t p);
void pst_img_graph_write_edge(FILE* arq, pst_img_graph_t *g, oct_arc_t e);
void pst_img_graph_write(FILE* arq, pst_img_graph_t* g);

void pst_img_graph_read_vertex(FILE* arq, pst_img_graph_t * g);
pst_path_t pst_img_graph_read_path(FILE* arq);
void pst_img_graph_read_edge(FILE* arq, pst_img_graph_t *g, long int list_onext[], int list_long_bit[]);
pst_img_graph_t* pst_img_graph_read(FILE* arq);
long int* pst_img_graph_sort_equations(pst_imgsys_t *S,long int* ref_tab,double *iW ,long int g_n );

double pst_img_graph_compute_left_face_curl(pst_img_graph_t* g, oct_arc_t e0);
r2_t pst_img_graph_left_face_baricenter(pst_img_graph_t* g, oct_arc_t e0);
void pst_img_graph_put_curl_into_image(pst_img_graph_t* g,float_image_t* OZ);

void debug_vertex_remove_edge(pst_img_graph_t* g, oct_arc_t e);

void pst_img_graph_vertex_remove_general(
  pst_img_graph_t* g,
  long int vi,
  long int n_neighbours,
  double* w_i,
  double wmag,
  bool_t merge_diagonals,
  bool_t verbose);
