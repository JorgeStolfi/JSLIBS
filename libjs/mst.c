/* See rn_classif_mst.h. */
/* Last edited on 2023-03-18 11:17:42 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>
#include <mst.h>

void mst_build_complete(int32_t n, mst_arc_cost_t *acost, int32_t P[], double C[], bool_t verbose)
  { 
    /* Brute force implementation of Prim's algorithm with bubble sorted queue. */
    if (verbose) { fprintf(stderr, "computing minimum spanning forest for %d nodes\n", n); }
        
    /* Allocate a queue and put all vertices in it with null predecessors: */
    int32_t *Q = notnull(malloc(n*sizeof(int32_t)), "no mem"); /* Sample index queue. */
    { int32_t i;  for (i = 0; i < n; i++) { Q[i] = i; P[i] = i; C[i] = +INF; } }

    auto void resortQ(int32_t k0, int32_t k1);
      /* Rearranges {Q[k0..k1-1]} so that {C[Q[k0]]} is minimum. */
    
    /* Prim's loop: */
    int32_t m = 0; /* Num vertices already settled: */
    while (m < n)
      { /* Vertices {Q[0..m-1]} have been settled. */
        /* Vertices {Q[m..n-1]} are unsettled. */
        /* For every vertex {u} in {Q[m..n-1]}, {C[u]} is the min arc cost from {u} to {Q[0..m-1]}. */
        /* Ensure that the cheapest is {Q[m]}: */
        /* !!! Should use a heap for {Q[m..n-1]} with root at {Q[n-1]}. !!! */
        resortQ(m,n);
        /* Get the next vetex to be settled: */
        int32_t u = Q[m]; 
        double Cu = C[u];
        if (verbose) { fprintf(stderr, "  selecting %5d cost = %8.6f\n", u, Cu); }
        m++;
        /* Update the cost of all remaining items in {Q}: */
        int32_t k;
        for (k = m; k < n; k++)
          { int32_t v = Q[k]; 
            if (verbose) { fprintf(stderr, "    checking  %5d cost = %8.6f", v, C[v]); }
            double Cvu = acost(v, u); /* Cost of arc {(v,u)}. */
            if (Cvu < C[v])
              { if (verbose) { fprintf(stderr, " --> %8.6f\n", Cvu); }
                C[v] = Cvu; P[v] = u;
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

