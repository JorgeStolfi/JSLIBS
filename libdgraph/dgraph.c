/* See {dgraph.h} */
/* Last edited on 2022-10-20 06:15:25 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <bool.h>
#include <affirm.h>
#include <fget.h>
#include <spmat.h>
#include <spmat_io.h>

#include <dgraph.h>

#define dgraph_trivial_elem (FALSE)
#define dgraph_elem_is_trivial(X) ((X)==FALSE)

void dgraph_elem_write(FILE *wr, bool_t *valP);
void dgraph_elem_read(FILE *rd, bool_t *valP);

spmat_impl(dgraph_t, dgraph, bool_t);

void dgraph_elem_write(FILE *wr, bool_t *valP) { fputc((*valP ? 'T' : 'F'), wr); }
void dgraph_elem_read(FILE *rd, bool_t *valP) { (*valP) = fget_bool(rd); }

spmat_io_impl(dgraph_t, dgraph, bool_t);

/* PROCEDURES SPECIFIC TO GRAPHS: */

dgraph_vertex_count_t *dgraph_degrees(dgraph_t *G, int32_t which)
  { dgraph_vertex_count_t nv = (which == 0 ? G->rows : G->cols);
    dgraph_edge_count_t ne = G->ents;
    dgraph_vertex_count_t *deg = notnull(malloc(nv*sizeof(dgraph_vertex_index_t)), "no mem");
    dgraph_vertex_index_t v;
    for (v = 0; v < nv; v++) { deg[v] = 0; }
    dgraph_edge_index_t e; 
    for (e = 0; e < ne; e++)
      { if (G->e[e].val) 
          { dgraph_vertex_index_t v = (which == 0 ? G->e[e].row : G->e[e].col); 
            assert(v < nv);
            deg[v]++; 
          }
      }
    return deg;
  }
  
dgraph_edge_count_t dgraph_add_undirected_edge
  ( dgraph_t *G, 
    dgraph_vertex_index_t u, 
    dgraph_vertex_index_t v, 
    dgraph_edge_count_t ne
  )
  {
    ne = dgraph_add_element(G, ne, u, v, TRUE);
    if (u != v) { ne = dgraph_add_element(G, ne, v, u, TRUE); }
    assert(ne <= G->ents);
    return ne;
  }
 
dgraph_vertex_index_t *dgraph_find_spanning_forest(dgraph_t *G)
  {
    assert(G->cols == G->rows);
    dgraph_vertex_count_t nv = G->cols;
    dgraph_edge_count_t ne = G->ents;

    /* The spanning forest, initially trivial: */
    dgraph_vertex_index_t *parent = notnull(malloc(nv*sizeof(dgraph_vertex_index_t)), "no mem");
    dgraph_vertex_index_t v;
    for (v = 0; v < nv; v++) { parent[v] = v; }

    /* Scan the edges of {G} and merge trees by the union-find algorithm: */
    dgraph_edge_index_t e; 
    for (e = 0; e < ne; e++)
      { if (G->e[e].val) 
          { dgraph_vertex_index_t r1 = dgraph_find_root(parent, G->e[e].row);  
            dgraph_vertex_index_t r2 = dgraph_find_root(parent, G->e[e].col);  
            if (r1 != r2) { parent[r1] = r2; }
          }
      }

    return parent;
  }

dgraph_vertex_index_t dgraph_find_root(dgraph_vertex_index_t parent[], dgraph_vertex_index_t u)
  {
    /* Find the root {r}: */
    dgraph_vertex_index_t r = u; 
    dgraph_vertex_index_t s = parent[r];
    while (s != r) { r = s; s = parent[s]; }
    /* Flatten the path: */
    dgraph_vertex_index_t v = u;
    while (v != r)
      { dgraph_vertex_index_t w = parent[v];
        parent[v] = r;
        v = w;
      }
    return r;
  }

dgraph_vertex_count_t *dgraph_count_components_by_size(dgraph_t *G)
  {
    assert(G->cols == G->rows);
    dgraph_vertex_count_t nv = G->cols;

    /* Find a maximally connected spanning forest of {G}: */
    dgraph_vertex_index_t *parent = dgraph_find_spanning_forest(G);

    /* For each root {v}, compute the size {tsize[v]} of its tree: */
    dgraph_vertex_count_t *tsize = notnull(malloc(nv*sizeof(dgraph_vertex_count_t)), "no mem");
    dgraph_vertex_index_t v;
    for (v = 0; v < nv; v++) { tsize[v] = 0; }
    for (v = 0; v < nv; v++) 
      { dgraph_vertex_index_t r = dgraph_find_root(parent, v); 
        assert(r < nv);
        tsize[r]++;
      }

    /* Count components by size: */
    dgraph_vertex_count_t *ct = notnull(malloc((nv+1)*sizeof(dgraph_vertex_count_t)), "no mem");
    dgraph_vertex_count_t size;
    for (size = 0; size <= nv; size++) { ct[size] = 0; }
    for (v = 0; v < nv; v++) 
      { if (parent[v] == v)
          { /* Vertex {v} is a forest root; get the vertex count {size} of its tree: */
            size = tsize[v];
            assert(size > 0);
            assert(size <= nv);
            ct[size]++;
          }
      }

    /* Cleanup and return: */
    free(tsize);
    free(parent);
    return ct;
  }

