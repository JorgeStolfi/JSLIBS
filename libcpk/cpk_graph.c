/* See cpk_graph.h */
/* Last edited on 2024-12-31 14:43:56 by stolfi */ 

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <affirm.h>
#include <bool.h>

#include <cpk_graph.h>

cpk_graph_t cpk_gather_neighbors(uint32_t nV, ui2_vec_t E)
  {
    uint32_t nE = E.ne;

    /* Allocate areas: */
    cpk_graph_t G;
    G.nV = nV;
    G.nE = nE;
    G.nbr = talloc(3*nE, uint32_t);
    G.deg = talloc(nV, uint32_t);
    G.fnb = talloc(nV+1, uint32_t);
    
    /* Compute {G.deg[u]} for all {u}: */
    for (int32_t u = 0; u< nV; u++) { G.deg[u] = 0; }
    for (uint32_t j = 0; j < nE; j++) 
      { ui2_t *Ej = &(E.e[j]);
        uint32_t u = ORG(*Ej), v = DST(*Ej); 
        /* Sanity check: */
        assert((u >= 0) && (u < nV));
        assert((v >= 0) && (v < nV)); 
        /* No self-loops: */
        assert(u != v);
        G.deg[u]++; G.deg[v]++;
      }
    
    /* Compute {G.fnb[u] = SUM{G.deg[0..u-1]}, for {u=0..nV}: */
    { uint32_t k = 0;
      for (int32_t u = 0; u < nV; u++) { G.fnb[u] = k; k += G.deg[u]; }
      G.fnb[nV]=k;
      assert(k==2*nE);
    }
 
    /* Gather neighbors of each vertex: */
    for (uint32_t u = 0; u< nV; u++) { G.deg[u] = 0; }
    for (uint32_t j = 0; j < nE; j++) 
      { ui2_t *Ej = &(E.e[j]);
        uint32_t u = ORG(*Ej), v = DST(*Ej); 
        assert(u != v);
        G.nbr[G.fnb[u]+G.deg[u]] = v; G.deg[u]++; 
        G.nbr[G.fnb[v]+G.deg[v]] = u; G.deg[v]++;
      }
    return G;
  }

bool_t cpk_graph_adjacent(cpk_graph_t *G, uint32_t u, uint32_t v)
  { if (u == v) { return FALSE; /* Graphs are assumed loopfree. */ } 
    uint32_t du = G->deg[u];
    uint32_t dv = G->deg[v];
    if (du < dv)
      { /* Look for {v} among the neighbors of {u}: */
        uint32_t *nbu = &(G->nbr[G->fnb[u]]); 
        for (int32_t i = 0; i < du; i++) { if (nbu[i] == v) { return TRUE; } } 
      }
    else
      { /* Look for {u} among the neighbors of {v}: */
        uint32_t *nbv = &(G->nbr[G->fnb[v]]);
        for (int32_t i = 0; i < dv; i++) { if (nbv[i] == u) { return TRUE; } } 
      }
    return FALSE;
  }

uint32_t cpk_graph_count_common_neighbors(cpk_graph_t *G, uint32_t u, uint32_t v, bool_t temp[])
  { uint32_t du = G->deg[u];
    uint32_t dv = G->deg[v];
    if (u == v) { return du; } 
    uint32_t count = 0;
    /* Mark all neighbors of {u}: */
    uint32_t *nbu = &(G->nbr[G->fnb[u]]); 
    for (int32_t i = 0; i < du; i++) { uint32_t w = nbu[i]; assert(!temp[w]); temp[w] = TRUE; } 
    /* Count neighbors of {v} that are marked: */
    uint32_t *nbv = &(G->nbr[G->fnb[v]]);
    for (int32_t i = 0; i < dv; i++) { if (temp[nbv[i]]) { count++; } } 
    /* Be a good boy and clean up after yourself: */
    for (int32_t i = 0; i < du; i++) { temp[nbu[i]] = FALSE; } 
    return count;
  }

