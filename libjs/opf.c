/* See rn_classif_opf.h. */
/* Last edited on 2023-03-18 11:16:48 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <opf.h>

void opf_build_complete(int32_t n, double C[], opf_arc_cost_t *acost, opf_path_cost_t *pcost, int32_t P[], int32_t R[], bool_t verbose)
  { 
    /* Implementation with brute-force min find. */
    if (verbose) { fprintf(stderr, "computing optimum path forest for %d nodes\n", n); }
        
    /* Allocate a queue and put all vertices in it with null predecessors: */
    int32_t *Q = notnull(malloc(n*sizeof(int32_t)), "no mem"); /* Sample index queue. */
    { int32_t u;  for (u = 0; u < n; u++) { Q[u] = u; P[u] = u; } }

    auto void resortQ(int32_t k0, int32_t k1);
      /* Rearranges {Q[k0..k1-1]} so that {C[Q[k0]]} is minimum. */
    
    /* Dijkstra's loop: */
    int32_t m = 0; /* Num vertices already settled: */
    while (m < n)
      { /* Vertices {Q[0..m-1]} have been settled and are sorted by increasing {C}. */
        /* Vertices {Q[m..n-1]} are unsettled. */
        /* Ensure that the cheapest is {Q[m]}: */
        /* !!! Should use a heap for {Q[m..n-1]} with root at {Q[n-1]}. !!! */
        resortQ(m,n);
        /* Get the next vetex to be settled: */
        int32_t u = Q[m]; 
        double Cup = C[u];
        if (verbose) { fprintf(stderr, "  selecting %5d cost = %8.6f\n", u, Cup); }
        assert(C[P[u]] <= C[u]);
        /* Propagate the root label: */
        R[u] = (P[u] == u ? u : R[P[u]]);
        m++;
        /* Update the cost of all remaining items in {Q}: */
        int32_t k;
        for (k = m; k < n; k++)
          { int32_t v = Q[k]; 
            if (verbose) { fprintf(stderr, "    checking  %5d cost = %8.6f", v, C[v]); }
            double Ca = acost(v, u); /* Cost of arc {(v,u)}. */
            double Cvup = pcost(Ca, Cup); /* Cost of any optimum path that starts with {(v,u,...)} */
            assert(Cvup >= Cup);
            if (Cvup < C[v])
              { if (verbose) { fprintf(stderr, " --> %8.6f\n", Cvup); }
                C[v] = Cvup; P[v] = u;
              }
            else
              { if (verbose) { fprintf(stderr, "\n"); } }
          }
      }
    free(Q);

    /* Local procedure implementations */
    
    void resortQ(int32_t k0, int32_t k1)
      { assert(k0 < k1);
        int32_t t = k0;
        int32_t u = Q[t];
        double Cup = C[u];
        int32_t k;
        for (k = k0 + 1; k < k1; k++) 
          { if (C[Q[k]] < Cup) { t = k; u = Q[t]; Cup = C[u]; } }
        if (t != k0) { Q[t] = Q[k0]; Q[k0] = u; } 
      }
  }

