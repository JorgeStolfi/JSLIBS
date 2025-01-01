/* See cpk_mis.h */
/* Last edited on 2024-12-31 16:28:14 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <pqueue.h>

#include <cpk_mis.h>
#include <cpk_build.h>
#include <cpk_graph.h>
#include <cpk_io.h>
#include <cpk_debug.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */
  
cpk_mis_state_t *cpk_mis_state_new(cpk_graph_t *G, double W[])
  {
    uint32_t nV = G->nV;
    cpk_mis_state_t *SQ;
    SQ = (cpk_mis_state_t *)notnull(malloc(sizeof(cpk_mis_state_t)), "no mem");
    /* Externally allocated components: */
    SQ->G = G;
    SQ->W = W;
    /* Internally allocated components: */
    SQ->S =    talloc(nV, uint32_t);
    SQ->degS = talloc(nV, uint32_t);
    SQ->locS = talloc(nV, uint32_t);
    SQ->Q = pqueue_new();
    pqueue_realloc(SQ->Q, /*nmax*/ nV, /*zlim:*/ nV);
    pqueue_set_order(SQ->Q, -1);
    SQ->degQ = talloc(nV, uint32_t);
    /* Initialize tables: */
    SQ->nS = 0;
    SQ->WS = 0;
    for (int32_t v = 0; v < nV; v++) { SQ->degS[v] = 0; SQ->locS[v] = locNONE; SQ->degQ[v] = 0; }
    return SQ;
  }

void cpk_mis_state_free(cpk_mis_state_t *SQ)
  {
    free(SQ->S);
    free(SQ->degS);
    free(SQ->locS);
    pqueue_free(SQ->Q);
    free(SQ->degQ);
    free(SQ);
  }

double cpk_mis_score(cpk_mis_state_t *SQ, uint32_t w)
  { 
    return SQ->W[w]/(double)(SQ->degQ[w] + 1);
  }

void cpk_mis_reset(cpk_mis_state_t *SQ, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "  Making S = {}, Q = {}\n"); }
    uint32_t nV = SQ->G->nV;
    SQ->nS = 0;
    SQ->WS = 0;
    pqueue_reset(SQ->Q);
    for (int32_t v = 0; v < nV; v++) 
      { SQ->degS[v] = 0;
        SQ->locS[v] = locNONE;
        SQ->degQ[v] = 0;
      }
  }

void cpk_mis_init(cpk_mis_state_t *SQ, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "  Making S = {}, Q = VG\n"); }
    uint32_t nV = SQ->G->nV;
    SQ->nS = 0;
    SQ->WS = 0;
    pqueue_reset(SQ->Q);
    for (uint32_t u = 0; u < nV; u++) 
      { SQ->degS[u] = 0;
        SQ->locS[u] = locNONE;
        SQ->degQ[u] = SQ->G->deg[u];
        double su = cpk_mis_score(SQ, u);
        pqueue_insert(SQ->Q, u, su);
      } 
  }

void cpk_mis_add_indep(cpk_mis_state_t *SQ, uint32_t u, bool_t verbose)
  {
    if (SQ->locS[u] == locNONE)
      { if (verbose) { fprintf(stderr, " S+%d", u); }
        /* Remove {u} from {Q}, if it is there: */
        cpk_mis_simply_delete_adm(SQ, u);
        /* Get neighbors of {u}: */
        uint32_t du = SQ->G->deg[u]; 
        uint32_t *nbu = &(SQ->G->nbr[SQ->G->fnb[u]]);
        /* Remove from {Q} any neighbor {w} of {u} that is there: */
        for (int32_t i = 0; i < du; i++) 
          { uint32_t w = nbu[i]; cpk_mis_simply_delete_adm(SQ, w); }
        /* Add {u} to {S}: */
        cpk_mis_simply_add_indep(SQ, u);
        assert(SQ->S[SQ->nS-1] == u);
      }
  }

void cpk_mis_delete_indep(cpk_mis_state_t *SQ, uint32_t u, bool_t verbose)
  {
    uint32_t k = SQ->locS[u];
    if (k != locNONE)
      { if (verbose) { fprintf(stderr, " S-%d", u); }
        assert(SQ->degS[u] == 0);  /* {u} should not be adjacent to {S}. */
        assert(! pqueue_has(SQ->Q, u));  /* {u} should not be in {Q}. */
        /* Remove {u} from {S}: */
        cpk_mis_simply_delete_indep(SQ, u);
        /* Add {u} to {Q}: */
        cpk_mis_simply_add_adm(SQ, u);
        /* Get neighbors of {u}: */
        uint32_t du = SQ->G->deg[u]; 
        uint32_t *nbu = &(SQ->G->nbr[SQ->G->fnb[u]]);
        /* Add to {Q} any neighbors that became admissible: */
        for (int32_t i = 0; i < du; i++)
          { uint32_t w = nbu[i];
            assert(SQ->locS[w] == locNONE); /* Vertex {w} must not be {S}. */
            /* If {w} became admissible, add it to {Q}: */
            if ((SQ->degS[w] == 0) && (! pqueue_has(SQ->Q, w)))
              { cpk_mis_simply_add_adm(SQ, w); }
          }
      }
  }

void cpk_mis_copy_indep(cpk_mis_state_t *SQ, uint32_t nX, uint32_t X[], bool_t verbose)
  {
    cpk_mis_init(SQ, verbose);
    for (int32_t i = 0; i < nX; i++) { cpk_mis_add_indep(SQ, X[i], verbose); }
  }
  
void cpk_mis_simply_add_adm(cpk_mis_state_t *SQ, uint32_t u)
  { 
    if (! pqueue_has(SQ->Q, u))
      { assert(SQ->locS[u] == locNONE);  /* Vertex {u} must not be in {S}. */
        assert(SQ->degS[u] == 0);   /* Vertex {u} must not be adjacent to {S}. */
        /* if (VERBOSE) { fprintf(stderr, " Q+%d", u); } */
        /* Add {u} to {Q}: */
        { double su = cpk_mis_score(SQ, u);
          pqueue_insert(SQ->Q, u, su);
        }
        /* Get neighbors of {u}: */
        uint32_t du = SQ->G->deg[u]; 
        uint32_t *nbu = &(SQ->G->nbr[SQ->G->fnb[u]]);
        /* Increment {degQ[w]} for every neighbor {w} of {u}: */
        for (int32_t i = 0; i < du; i++)
          { uint32_t w = nbu[i];
            assert(SQ->locS[w] == locNONE);  /* Vertex {w} must not be in {S}. */
            SQ->degQ[w]++;
            if (pqueue_has(SQ->Q, w))
              { /* Vertex {w} is in {Q}, recompute its score in the queue: */
                assert(SQ->degS[w] == 0);  /* Vertex {w} must not be adjacent to {S}. */
                /* Weight we would gain by taking {w}, divided by adms we would lose. */
                double sw = cpk_mis_score(SQ, w); 
                pqueue_set_value(SQ->Q, w, sw);
              }
          }
      }
  }
 
void cpk_mis_simply_add_indep(cpk_mis_state_t *SQ, uint32_t u)
  {
    if (SQ->locS[u] == locNONE)
      { assert(! pqueue_has(SQ->Q, u));  /* Vertex {u} should not be in {Q}. */
        assert(SQ->degQ[u] == 0);  /* Vertex {u} should not be adjacent to {Q}. */
        /* Get neighbors of {u}: */
        uint32_t du = SQ->G->deg[u]; 
        uint32_t *nbu = &(SQ->G->nbr[SQ->G->fnb[u]]);
        /* Add {u} to {S}: */
        assert(SQ->nS >= 0);
        SQ->S[SQ->nS] = u; SQ->locS[u] = SQ->nS; SQ->nS++; 
        /* Update weight: */
        SQ->WS += SQ->W[u];
        /* Update {degS[w]} for every neighbor {w} of {u}: */
        for (int32_t i = 0; i < du; i++) 
          { uint32_t w = nbu[i]; 
            assert(! pqueue_has(SQ->Q, w)); /* Vertex {w} must not be in {Q}. */
            SQ->degS[w]++;
          }
        assert(SQ->S[SQ->nS-1] == u);
      }
  }

void cpk_mis_simply_delete_indep(cpk_mis_state_t *SQ, uint32_t u)
  {
    uint32_t k = SQ->locS[u];
    if (k != locNONE)
      { assert(! pqueue_has(SQ->Q, u));  /* Vertex {u} should not be in {Q}. */
        assert(SQ->degQ[u] == 0);  /* Vertex {u} should not be adjacent to {Q}. */
        /* Remove {u} from {S}: */
        { SQ->nS--; uint32_t v = SQ->S[SQ->nS]; SQ->S[k] = v; SQ->locS[v] = k; SQ->locS[u] = locNONE; }
        /* Update weight: */
        SQ->WS -= SQ->W[u];
        /* Get neighbors of {u}: */
        uint32_t du = SQ->G->deg[u]; 
        uint32_t *nbu = &(SQ->G->nbr[SQ->G->fnb[u]]);
        /* Decrement {degS[w]} for every neighbor {w} of {u}: */
        for (int32_t i = 0; i < du; i++)
          { uint32_t w = nbu[i];
            assert(SQ->locS[w] == locNONE); /* {w} should not be in {S}: */
            assert(! pqueue_has(SQ->Q, w)); /* {w} should not be in {Q}: */
            assert(SQ->degS[w] > 0); /* {w} was adjacent to a vertex of {S}: */
            /* decrement {S} degree: */
            SQ->degS[w]--;
          }
      }
  }

void cpk_mis_simply_delete_adm(cpk_mis_state_t *SQ, uint32_t u)
  { if (pqueue_has(SQ->Q, u)) 
      { /* fprintf(stderr, " Q-%d", u); */
        /* Consistency checks: */
        assert(SQ->locS[u] == locNONE);
        assert(SQ->degS[u] == 0);
        /* if (VERBOSE) { fprintf(stderr, " Q-%d", u); } */
        /* Remove {u} from {Q} */
        pqueue_delete(SQ->Q, u);
        /* Decrement {degQ[w]} for every neighbor {w} of {u}: */
        uint32_t du = SQ->G->deg[u]; 
        uint32_t *nbu = &(SQ->G->nbr[SQ->G->fnb[u]]);
        for (int32_t i = 0; i < du; i++)
          { uint32_t w = nbu[i];
            SQ->degQ[w]--; assert(SQ->degQ[w] >= 0);
            /* Recompute the score of {w} in the queue: */
            if (pqueue_has(SQ->Q, w))
              { /* Consistency checks: */
                assert(SQ->locS[w] == locNONE);
                assert(SQ->degS[w] == 0);
                /* Weight we would gain by taking {w}, divided by adms we would lose. */
                double sw = cpk_mis_score(SQ, w); 
                pqueue_set_value(SQ->Q, w, sw);
              }
          }
      }
  }
