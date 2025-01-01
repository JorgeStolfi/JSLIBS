/* See cpk_greedy.h */
/* Last edited on 2024-12-31 14:45:18 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <bool.h>
#include <pqueue.h>
#include <affirm.h>
#include <jsrandom.h>

#include <cpk_greedy.h>
#include <cpk_graph.h>
#include <cpk_mis.h>
#include <cpk_lopt.h>
#include <cpk_io.h>
#include <cpk_debug.h>

/* INTERNAL PROTOTYPES */

uint32_t cpk_count_common_Q_neighbors(cpk_graph_t *G, pqueue_t *Q, uint32_t u, uint32_t v, bool_t temp[]);
  /* Returns the number of vertices of {Q} that are adjacent to both {u} and {v}. */

/* IMPLEMENTATIONS */

uint32_vec_t cpk_greedy_find_indep_set
  ( cpk_graph_t *G,  /* The incompatibilty graph. */
    double W[],      /* Weight of each vertex. */
    bool_t verbose   /* TRUE prints debugging diagnostics. */
  )
  {
    if (verbose) { cpk_trace_entry(); }

    /* Call generic greedy search procedure: */
    cpk_mis_state_t *SQ = cpk_mis_state_new(G, W);
    cpk_mis_init(SQ, verbose);

    if (verbose) { fprintf(stderr, "\n------------------------------------\n"); }
    if (verbose) { fprintf(stderr, "Greedily finding an independent set"); }
    cpk_greedy_maximize(SQ, verbose);
    if (verbose) { fprintf(stderr, "\n"); }
    if (verbose) { TRACE_SOL("Greedy solution", SQ->nS, SQ->WS);  }
   
    /* Save it as the seen so far (for local opt): */
    uint32_t SM[G->nV]; 
    uint32_t nSM = 0;
    double WSM = 0;
    cpk_lopt_update_best_solution(SQ, &nSM, SM, &WSM, verbose);
    
    if (verbose) { fprintf(stderr, "\n------------------------------------\n"); }
    if (verbose) { fprintf(stderr, "Applying local optimization"); }
    uint32_t maxRem = 2;
    cpk_local_opt(SQ, maxRem, &nSM, SM, &WSM, verbose);
    if (verbose) { fprintf(stderr, "\n"); }
    
    /* Copy the best solution into the current one: */
    cpk_mis_copy_indep(SQ, nSM, SM, verbose);
    affirm(cpk_weight_cmp(SQ->WS, WSM) == 0, "bug");
    if (verbose) { TRACE_SOL("Locally optimized solution", SQ->nS, SQ->WS);  }
    
    /* Repackage the independent set for return: */
    uint32_vec_t J = uint32_vec_new(nSM);
    for (int32_t i = 0; i < nSM; i++) { J.e[i] = SM[i]; }
    
    /* Housecleaning: */
    cpk_mis_state_free(SQ);
    
    if (verbose) { cpk_trace_exit(); }
    return J;
  }

void cpk_greedy_maximize(cpk_mis_state_t *SQ, bool_t verbose)
  { 
    while (pqueue_count(SQ->Q) > 0)
      { uint32_t u = pqueue_head(SQ->Q); 
        cpk_mis_add_indep(SQ, u, verbose);
      }
  }

void cpk_greedy_maximize_random(cpk_mis_state_t *SQ, bool_t verbose)
  { uint32_t nQ;
    pqueue_t *Q = SQ->Q;
    while ((nQ = pqueue_count(Q)) > 0)
      { /* Pick a random index in level {r} of the heap with probability 3*(1/4)^r: */
        uint32_t a = 0, b = 1;
        double r = drandom();
        while ((r < 0.25) && (b < nQ)) { a = b; b = 2*a+1; r *= 4; }
        if (b > nQ) { b = nQ; }
        uint32_t i = uint32_abrandom(a, b-1);
        if (verbose) { fprintf(stderr, " [%d]", i); }
        pqueue_item_t u = pqueue_item(Q, i); 
        cpk_mis_add_indep(SQ, u, verbose);
      }
    if (verbose) { fprintf(stderr, "\n");  }
  }

void cpk_greedy_maximize_lookahead(cpk_mis_state_t *SQ, bool_t verbose, bool_t temp[])
  { cpk_graph_t *G = SQ->G;
    pqueue_t *Q = SQ->Q;
    uint32_t nQ;
    while ((nQ = pqueue_count(Q)) > 0)
      { /* procura par independente {u,v} em {Q[0..nQ-1]} com melhor escore conjunto: */
        double suv = -INF; 
        uint32_t u = UINT32_MAX; uint32_t v = UINT32_MAX;
        for (uint32_t ix = 0; ix < nQ; ix++)
          { pqueue_item_t x = pqueue_item(Q, ix);
            for (uint32_t iy = ix+1; iy < nQ; iy++)
              { pqueue_item_t y = pqueue_item(Q, iy);
                if (! cpk_graph_adjacent(G, x, y))
                  { /* Conta vizinhos comuns: */
                    uint32_t dxyQ = cpk_count_common_Q_neighbors(G, Q, x, y, temp);
                    /* Calcula razão (peso ganho)/(candidatos perdidos) do par: */
                    double sxy = (double)(SQ->W[x] + SQ->W[y])/((double)(SQ->degQ[x] + SQ->degQ[y] + 2) - dxyQ); 
                    /* Guarda melhor par: */
                    if (sxy > suv) { suv = sxy; u = x; v = y; }
                  }
              }
          }
        assert((u < nQ) && (v < nQ));
        
        /* Escolhe melhor {w} para incluir em {S}: */
        uint32_t w;
        if (suv == -INF)
          { /* Não há mais pares independentes em {Q}. */
            w = pqueue_head(Q);
          }
        else
          { /* Escolhe {w} como o melhor dentre {u,v}: */
            double su = SQ->W[u]/(G->deg[u] + 1);
            double sv = SQ->W[v]/(G->deg[v] + 1);
            w = (su >= sv ? u : v);
          }
        
        /* Acrescenta {w} a {S}, e atualiza peso: */
        cpk_mis_add_indep(SQ, w, verbose);
      }
  }

uint32_t cpk_count_common_Q_neighbors(cpk_graph_t *G, pqueue_t *Q, uint32_t u, uint32_t v, bool_t temp[])
  { uint32_t du = G->deg[u];
    uint32_t dv = G->deg[v];
    if (u == v) { return du; } 
    uint32_t count = 0;
    /* Mark all neighbors of {u}: */
    uint32_t *nbu = &(G->nbr[G->fnb[u]]); 
    for (int32_t i = 0; i < du; i++) 
      { uint32_t w = nbu[i]; assert(!temp[w]); temp[w] = TRUE; } 
    /* Count neighbors of {v} that are marked and in {Q}: */
    uint32_t *nbv = &(G->nbr[G->fnb[v]]);
    for (int32_t i = 0; i < dv; i++) 
      { uint32_t w = nbv[i]; if (temp[w] && pqueue_has(Q,w)) { count++; } } 
    /* Be a good boy and clean up after yourself: */
    for (int32_t i = 0; i < du; i++) 
      { uint32_t w = nbu[i]; assert(temp[w]); temp[w] = FALSE; } 
    return count;
  }
