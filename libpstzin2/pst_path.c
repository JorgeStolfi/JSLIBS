/* See {pst_path.h} */
/* Last edited on 2024-12-23 13:53:32 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <r2.h>

#include <pst_path.h>

pst_path_t pst_path_create_empty(void)
  { pst_path_t p = (pst_path_t){ .n=0, .v=NULL, .reverse=FALSE }; 
    return p;
  }
  
void pst_path_free(pst_path_t p)
  { free(p.v); }

pst_path_t pst_path_create_single(r2_t coords)
  { r2_t *v = (r2_t*)malloc(sizeof(r2_t));
    v[0] = coords;
    pst_path_t p = (pst_path_t){ .n=1, .v=v, .reverse=FALSE }; 
    return p;
  }

pst_path_t pst_path_reverse(pst_path_t p)
  { p.reverse = !p.reverse;
    return p;
  }

r2_t pst_path_get_vertex(pst_path_t p, uint32_t i)
  { return ( p.reverse ? p.v[p.n - i -1 ] : p.v[i] );  }

pst_path_t pst_path_concatenate(pst_path_t p0, r2_t coords, pst_path_t p1)
  { uint32_t n0 = p0.n;
    uint32_t n1 = p1.n;
    uint32_t n = n0 + n1 + 1;

    r2_t *v = (r2_t*)malloc(sizeof(r2_t)*n);
    for (uint32_t i = 0; i < n; i++)
      { r2_t c;
        if (i < n0)
          { c = pst_path_get_vertex(p0,i); }
        else if (i == n0)
          { c = coords; }
        else
          { uint32_t j = (uint32_t)(i - n0 - 1);
            c = pst_path_get_vertex(p1,j);
          }
        v[i] = c;
      }
    return (pst_path_t) { .n=n, .v=v, .reverse=FALSE };
  }

void pst_img_graph_write_path(FILE *wr, pst_path_t p)
  {
    fprintf(wr,"%d %d ",p.n,p.reverse);
    for (uint32_t i = 0; i < p.n; i++)
      { fprintf(wr,"%9.6f %9.6f ",p.v[i].c[0],p.v[i].c[1]); }
  }

pst_path_t pst_img_graph_read_path(FILE *wr)
  {
    pst_path_t p;
    int32_t reverse;
    demand(fscanf(wr,"%d %d",&(p.n),&(reverse)) == 2, "Cannot read path");
    p.reverse = (reverse == 1);
    p.v = NULL;
    if (p.n == 0) return p;
    p.v = (r2_t*) malloc(sizeof(r2_t)*(p.n));
    for (uint32_t i = 0; i < p.n; i++)
      { int32_t nr1 = fscanf(wr,"%lf %lf",&(p.v[i].c[0]),&(p.v[i].c[1]));
        demand(nr1 == 2, "Cannot read path elements");
      }
    return p;
  }
