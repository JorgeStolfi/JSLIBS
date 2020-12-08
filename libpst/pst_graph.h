#ifndef pst_graph_H
#define pst_graph_H

/* Last edited on 2013-03-18 01:40:20 by stolfilocal */
/* Created by Rafael F. V. Saracchini */

#include <float_image.h>
#include <pst_imgsys.h>


typedef struct pst_edge_t{
  /* Represents an edge u->v */
  long int u,v; /* It stores the index on the array of vertices (u,v) */
  double d; /* Derivative d of the edge */
  double w; /*Weight w of the edge */
  int axis; /*Used to identify the orientation of the edge, notice that altought meaningless in the 
	      projection plane, it is usefull to determine 4-neighbours of a vertex*/
} pst_edge_t;

typedef struct pst_vertex_t{
  /* Represents a vertex */
  long int id; /* Stores unique id number for each vertex, to allow identification */
  long int i,j; /* Stores vertex index, used only for removal purposes !*/
  
  long int vtxm,vtxp,vtym,vtyp; /* Stores the edge's indexes that connects the 4 neighbours of the vertex. It will contain -1 if the
				   edge does not exist*/
} pst_vertex_t;

struct pst_graph_t{
  pst_vertex_t* vertices; /* Array that contain the vertices of the graph*/
  pst_edge_t* edges; /* Array that contain the edges of the graph*/
  long int num_vertex; /* How much vertices exists in g*/ 
  long int num_edge; /* How much edges exists in g */
  long int max_vertex; /* Maximum allowed number of vertices*/
  long int max_edge; /* Maximum number allowed of edges */
};

typedef struct pst_graph_t pst_graph_t;

pst_graph_t* pst_graph_new(long int max_vertex, long int max_edge);
/*Alocate graph structure with enough space for {max_vertex} vertices and {max_edge} edges*/


void pst_graph_add_vertex(pst_graph_t* g, long int id, long int i, long int j);
/* Add an vertex with id {id}, and indexes {(i,j)}. Initializes the vertex as it was completely disconected from the graph*/

void pst_graph_add_edge(pst_graph_t* g, long int u, long int v, double d, double w, int axis);
/* Add an edge to the graph from {u} to {v} with  derivative {d} and weight {w}. It allows add a edge to vertices
that arent still added to the graph, and does not update the neighbours of the vertices because this.  */

void pst_graph_update_neighbours(pst_graph_t* g);
/* Updates the neighbours of all vertices that have an edge in the graph, it will fail if the edges arent consistent wuth the
vertices.*/



long int pst_graph_compute_vertex_index(long int ix, long int iy, long int NX, long int NY);
/* Given a index{(ix,iy)} it computes a unique id, given the dimensions {NX}x{NY} of the original image*/

void pst_graph_restore_vertex_index(long int id, long int NX, long int NY, long int *ix, long int *iy);
/* Returns a {(ix,iy)} index from a unique id computed by pst_graph_compute_vertex_index given the original image dimensions {NX}x{NY} */

void pst_graph_get_edge_derivate_weight(pst_graph_t* g,long int v, long int e,double*d, double* w);
/* Given an vertex{v} and a edge{e} returns the equivalent derivate {d} and weight {w}. If {v} is the origin of the
edge, it returns e.d, if it is the endpoint, returns -e.d. It fails if v does not belong to the edge*/

void pst_graph_vertex_get_neighbour(pst_graph_t *g,
				    long int vt0,
				    int dx, int dy,
				    long int *vt1,
				    double *d,
				    double *w
				    );
/* Given a vertex {vt0} and a direction {(dx,dy)} where -1 <= dx <=1  and -1<=dy<=1 returns a vertex {vt1} in that direction and
the derivate {d} and weight{w} of the path. If there isnt a simple path to such edge, it returns -1 as value of {vt1} and {d} and
{w} values will be 0. If fails if the given direction is bigger than allowed or if it is (0,0) */
				    


pst_graph_t* pst_graph_create_from_gradient_weight_map(float_image_t* IG, float_image_t* IW);
/* Create a graph from a gradient image {IG} and weight image ${IW} */


void pst_graph_write_vertex(FILE* arq, pst_vertex_t* v);
/* Write in the stream {arq} the vertex {v} in human-readable fashion*/
void pst_graph_write_edge(FILE* arq, pst_edge_t* e);
/* Write in the stream {arq} the edge{e} in human-readable fashion*/
void pst_graph_write(FILE* arq,pst_graph_t* g);
/* Write in the stream {arq} the graph {g}  using the pst_graph_write_vertex and pst_graph_write_edge*/
bool_t pst_graph_check_consistency(FILE* arq,pst_graph_t* g);
/* Check the graph consistensy using the edge's list. If {arq} is not null, a relatory is written in humam-readabla fashin   */

pst_graph_t* pst_graph_shrink(pst_graph_t* g);
/* Shrinks the graph...*/

int* pst_graph_compute_mst(pst_graph_t* g);

void pst_graph_interpolate_two_samples
  (  float_image_t* I, float_image_t* W,
     int c,
     int x0, int y0,
     int x1, int y1,
     double *v, double* w
   );

void pst_graph_mst_walk_recursive(pst_graph_t* g, long int vertex,int* marked_vertices, int* marked_edges,int num_tree);

pst_imgsys_t* pst_graph_build_integration_system(pst_graph_t* g,long int NX_Z,long int NY_Z);

void pst_graph_solve_system(
    pst_graph_t* g,
    pst_imgsys_t* S,
    float_image_t *OZ, 
    long int maxIter, 
    double convTol, 
    int para, 
    int szero, 
    bool_t verbose 
  );
  
void pst_graph_estimate_from_shrunk(pst_graph_t* g, pst_graph_t* jg, float_image_t* OZ);

void pst_graph_integration_recursive( 
  pst_graph_t* g,
  float_image_t* OZ,
  long int maxIter,
  double convTol, 
  int para, 
  int szero, 
  bool_t verbose,
  int level
);

void pst_graph_free(pst_graph_t* g);

#endif
