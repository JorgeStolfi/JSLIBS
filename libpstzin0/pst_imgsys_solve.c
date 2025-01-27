/* See pst_imgsys_solve.h  */

/* Created on 2005-10-01 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-01-25 08:46:08 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <rn.h>
#include <affirm.h>
#include <jsmath.h>

#include <pst_imgsys.h>

#include <pst_imgsys_solve.h>

double pst_imgsys_sol_change(double *z_old, double *z, uint32_t N);
  /* Returns the maximum of {abs(z_old[k]-z[k])}, for {k = 0..N-1}. */

void pst_imgsys_solve_iterative
  ( pst_imgsys_t *S, 
    double z[],
    uint32_t ord[],
    uint32_t maxIter, 
    double convTol,
    bool_t para, 
    bool_t szero,
    bool_t verbose,
    int32_t level,
    uint32_t reportStep,
    pst_imgsys_solve_report_sol_proc_t *reportSol
  )
  {
    uint32_t N = S->N;
    int32_t indent = (level < -1 ? 0 : 2*level+2);

    /* Previous solution: */
    double *z_old = rn_alloc(N);
    
    double change = +INF;  /* Max {z} change in last iteration. */
    int32_t iter = 0;
    while (iter < maxIter)
      {
        /* Report the current solution if so requested: */
        if (reportSol != NULL)
          { bool_t reportSolBefore = (iter == 0);
            if (reportStep != 0)
              { reportSolBefore |= (iter < reportStep);
                reportSolBefore |= ((iter % (int32_t)reportStep) == 0);
              }
            if (reportSolBefore) { reportSol(level, iter, change, FALSE, N, z); }
          }

        /* Another pass over all variables {z[k]}: */
        for (uint32_t kk = 0; kk < N; kk++)
          { /* Choose the next variable to recompute: */
            uint32_t k = (ord != NULL ? ord[kk] : kk);
            /* Save current solution in {z_old}: */
            z_old[k] = z[k];
            /* Get equation {k}: */
            pst_imgsys_equation_t *eqk = &(S->eq[k]);
            /* Compute {z[k]} using equation {k}: */
            double sum = eqk->rhs; /* Right-hand side of equation. */
            double cf_k = 0.0; /* Coefficient of {z[k]} in equation. */
            for(uint32_t t = 0; t < eqk->nt; t++)
              { /* Get hold of another variable {z[j] entering in equation {k}: */
                uint32_t j = eqk->uid[t];
                demand((j >= 0) && (j < N), "invalid variable index in system");
                double cf_j = eqk->cf[t];
                if (j == k) 
                  { /* This term uses {z[k]}, store the coefficient: */
                    cf_k = cf_j;
                  }
                else if (cf_j != 0)
                  { /* The variable {z[j]} is distinct from {z[k]}. */
                    /* Get the appropriate value (new or old) of {z[j]}: */
                    double Zj = (para && (j < k) ? z_old[j] : z[j]);
                    /* Subtract from the right-hand side: */
                    sum = sum - Zj * cf_j;
                  }
              }
            
            /* Require that {eqk} depends on {z[k]}: */
            demand(cf_k != 0.0, "system's matrix has a zero in the diagonal"); 
            
            /* Solve the equation: */
            z[k] = sum / cf_k;
          }
          
        if (szero)
          { /* Normalize for zero sum: */
            double sum = 0;
            for (uint32_t k = 0; k < N; k++) { sum += z[k]; }
            double avg = sum/N;
            for (uint32_t k = 0; k < N; k++) { z[k] -= avg; }
          }
        iter++;
        
        /* Check for apparent convergence: */
        change = pst_imgsys_sol_change(z_old, z, N);
        if (verbose)
          { fprintf(stderr, "%*s iteration %3d change = %16.8f\n", indent, "", iter, change); }
        if (change <= convTol) { /* Converged: */ break; }
      }
      
    if (change > convTol)
      { /* Failed to converge: */
        fprintf(stderr, "%*s** gave up after %6d iterations, last change = %16.8f\n", indent, "", iter, change);
      }
    else
      { /* Converged: */
        fprintf(stderr, "%*sconverged after %6d iterations,last change = %16.8f\n", indent, "", iter, change);
      }

    if (reportSol != NULL) { reportSol(level, iter, change, TRUE, N, z); }

    free(z_old);
  }

double pst_imgsys_sol_change(double* z_old, double* z, uint32_t N)
  { uint32_t k;
    double max_change = 0;
    for(k = 0; k< N; k++)
      { double dk = z[k] - z_old[k];
        if(dk < 0) { dk = - dk; }
        if(dk > max_change) { max_change = dk; }
      }
    return max_change;
  }

uint32_t* pst_imgsys_sort_equations(pst_imgsys_t *S)
  {
    uint32_t M = 2*(pst_imgsys_MAX_COEFFS-1); /* Max number of arcs out of or into a node. */
    uint32_t N = S->N; /* Variables and equations are numbered from {0 to N-1}. */
    
    /* BUILDING THE PRIORITY GRAPH */

    /* The priority graph has one node per system variable,
      and at most one arc per non-diagonal coefficient of the system, 
      pointing from the variable with higher priority to that
      of lower priority. It may have duplicate arcs. */
    
    /* Arcs of the priority graph: */
    uint32_t *dst = talloc(M*N, uint32_t);
    /* Value of {n_out[k]} is count of neighbours of node {k} with lower priority. */
    uint32_t *n_out = talloc(N, uint32_t);
    /* Value of {n_in[k]} is count of unprocessed neighbours of node {k} with higher priority. */
    uint32_t *n_in = talloc(N, uint32_t);
    /* The arcs out of node {k} go to nodes {dst[M*k+j]} for {j} in {0..n_out[k]-1}. */

    auto void maybe_add_arrow(uint32_t k1, uint32_t k2);
      /* Adds an arc from node {k1} to node {k2}, 
        if the weight of {k1} is smaller than that of {k2}. */
    
    void maybe_add_arrow(uint32_t k1, uint32_t k2)
      { /* Nodes with lower weights have higher priority: */
        double w1 = S->eq[k1].wtot;
        double w2 = S->eq[k2].wtot;
        if (w1 < w2) { dst[M*k1 + n_out[k1]] = k2; n_out[k1]++;  n_in[k2]++; }
      }

    /* Gather the graph arcs: */
    uint32_t k;
    for(k = 0; k < N; k++){ n_in[k] = 0; n_out[k] = 0; }
    for(k = 0; k < N; k++)
      { pst_imgsys_equation_t *eqk = &(S->eq[k]);
        assert( eqk->uid[0] == k);
        uint32_t j;
        for(j = 1; j < eqk->nt; j++)
          { uint32_t i = eqk->uid[j];
            maybe_add_arrow(k,i);
            maybe_add_arrow(i,k);
          }
      }

    /* TOPOLOGICAL SORT OF THE PRIORITY GRAPH */
    
    /* A topological sort of the graph gives a list {ord[0..N-1]} 
      of nodes such that if {k1} has higher priority than {k2}
      then {k1} appears before {k2} in {ord}. */

    /* Queue of source nodes, and output ordering: */
    uint32_t *ord = talloc(N, uint32_t); /* Nodes in toporder. */
    uint32_t q_free = 0;  /* Index of first free entry in {ord}. */
    uint32_t q_start = 0; /* Index in {ord} of first unprocessed entry. */ 
    /* Nodes {ord[0..q_start-1]} have been sorted and have no outgoing or incoming edges. */
    /* Nodes {ord[q_start..q_free-1]} have been sorted and have no incoming edges. */

    auto void queue_insert(uint32_t index);
    void queue_insert(uint32_t index)
      { assert(q_free < N);
        ord[q_free] = index;
        q_free++;
      }

    auto uint32_t queue_remove(void);
    uint32_t queue_remove(void)
      { assert(q_start < q_free );
        uint32_t i = ord[q_start];
        q_start++;
        return i;
      }  

    /* Insert the nodes of highest priority: */
    for(k = 0; k < N; k++) { if (n_in[k] == 0) { queue_insert(k); } }

    while(q_start < q_free )
      { /* Get the next unprocessed node {k}: */
        uint32_t k = queue_remove();
        /* Add all its dominated neighbors to the queue, erase the arcs: */
        uint32_t j;
        for(j = 0; j < n_out[k]; j++)
          { /* Get a dominated neighbor {i}: */
            uint32_t i = dst[M*k + j];
            /* Erase the arc from {k} to {i}: */
            assert(n_in[i] > 0);
            n_in[i]--;
            /* If {i} becomes a source, add it to the queue: */
            if(n_in[i] == 0){ queue_insert(i); }
          }
        n_out[k] = 0;
      }
    assert(q_free == N);

    free(dst);
    free(n_in);
    free(n_out);

    return ord;
  }

/*
**
** Copyright (C) Jorge Stolfi, Unicamp.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty. Neither the author nor its employers are liable to
** any damages which may result from its use.
*/
