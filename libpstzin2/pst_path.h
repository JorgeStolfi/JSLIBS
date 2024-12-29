/* Last edited on 2024-12-23 14:51:47 by stolfi */
/* Created by Rafael F. V. Saracchini */

#ifndef pst_path_H
#define pst_path_H

/* Paths for drawing curved graph edges. */ 

#include <r2.h>

typedef struct pst_path_t
  { uint32_t n;     /* number of internal vertices*/
    r2_t* v;        /*coordinates of internal vertices*/
    bool_t reverse; /* TRUE if the vertices are stored in reversed manner*/
  } pst_path_t;

#define DEFAULT_WMAG 1.0
pst_path_t pst_path_create_empty(void);

void pst_path_free(pst_path_t p);

pst_path_t pst_path_create_single(r2_t coords);

pst_path_t pst_path_reverse(pst_path_t p);

r2_t pst_path_get_vertex(pst_path_t p, uint32_t i);

pst_path_t pst_path_concatenate(pst_path_t p0, r2_t coords, pst_path_t p1);
void pst_img_graph_write_path(FILE* wr, pst_path_t p);

pst_path_t pst_img_graph_read_path(FILE* rd);

#endif
