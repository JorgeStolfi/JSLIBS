/* See pst_imgsys.h  */

/* Created on 2005-10-01 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2013-05-24 04:29:49 by stolfilocal */
/* See the copyright and authorship notice at the end of this file.  */


#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include <pst_imgsys.h>
#include <pst_height_map.h>

#include <float_image.h>
#include <float_image_mscale.h>
#include <affirm.h>

double pst_imgsys_sol_change(double *Zold, double *Z, int N);
  /* Returns the maximum of {abs(Zold[k]-Z[k])}, for {k = 0..N-1}. */

pst_imgsys_t *pst_imgsys_new(int NX, int NY, int N, int *ix, int *col, int *row)
  {
    pst_imgsys_t *S = (pst_imgsys_t*)malloc(sizeof(pst_imgsys_t));
    S->N = N;
    S->eq = (pst_imgsys_equation_t*)malloc(sizeof(pst_imgsys_equation_t)*N);
    S->ix = ix;
    S->col = col;
    S->row = row;
    return S;
  }

pst_imgsys_t *pst_imgsys_from_eqs(int NX, int NY, int N, pst_imgsys_equation_t *eq, int *ix, int *col, int *row)
  {
    pst_imgsys_t *S = (pst_imgsys_t*)malloc(sizeof(pst_imgsys_t));
    S->N = N;
    S->eq = eq;
    S->ix = ix;
    S->col = col;
    S->row = row;
    return S;
  }

void pst_imgsys_free(pst_imgsys_t *S)
  {
    free(S->eq);
    free(S->col);
    free(S->row);
    free(S->ix);
    free(S);
  }

void pst_imgsys_solve
  ( pst_imgsys_t *S, 
    double Z[],
    int ord[],
    int maxIter, 
    double convTol,
    int para, 
    int szero,
    bool_t verbose,
    int indent,
    pst_imgsys_solution_report_proc_t *reportSol
  )
  {
    int N = S->N;

    /* Previous solution: */
    double *Zold = (double*)malloc(sizeof(double)*N);
    
    int iter;
    double change = +INF;  /* Max {Z} change in last iteration. */
    for (iter = 0; iter < maxIter; iter++)
      {
        /* Report the current solution if so requested: */
        if (reportSol != NULL) { reportSol(iter, change, FALSE, N, Z); }

        /* Another pass over all unknowns {Z[k]}: */
        int kk;
        for (kk = 0; kk < N; kk++)
          { /* Choose the next unknown to recompute: */
            int k = (ord != NULL ? ord[kk] : kk);
            /* Save current solution in {Zold}: */
            Zold[k] = Z[k];
            /* Get equation {k}: */
            pst_imgsys_equation_t *eqk = &(S->eq[k]);
            /* Compute {Z[k]} using equation {k}: */
            double sum = eqk->rhs; /* Right-hand side of equation. */
            double cf_k = 0.0; /* Coefficient of {Z[k]} in equation. */
            int t;
            for(t = 0; t < eqk->nt; t++)
              {
                /* Get hold of another unknown {Z[j] entering in equation {k}: */
                int j = eqk->ix[t];
                demand((j >= 0) && (j < N), "invalid variable index in system");
                double cf_j = eqk->cf[t];
                if (j == k) 
                  { /* This term uses {Z[k]}, store the coefficient: */
                    cf_k = cf_j;
                  }
                else if (cf_j != 0)
                  { /* The unknown {Z[j]} is distinct from {Z[k]}. */
                    /* Get the appropriate value (new or old) of {Z[j]}: */
                    double Zj = (para && (j < k) ? Zold[j] : Z[j]);
                    /* Subtract from the right-hand side: */
                    sum = sum - Zj * cf_j;
                  }
              }
            
            /* Require that {eqk} depends on {Z[k]}: */
            demand(cf_k != 0.0, "system's matrix has a zero in the diagonal"); 
            
            /* Solve the equation: */
            Z[k] = sum / cf_k;
          }
          
        if (szero)
          { /* Normalize for zero sum: */
            double sum = 0;
            int  k;
            for (k = 0; k < N; k++) { sum += Z[k]; }
            double avg = sum/N;
            for (k = 0; k < N; k++) { Z[k] -= avg; }
          }
        iter++;
        
        /* Check for apparent convergence: */
        change = pst_imgsys_sol_change(Zold, Z, N);
        if (verbose)
          { fprintf(stderr, "%*s  iteration %3d change = %16.8f\n", indent, "", iter, change); }
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

    if (reportSol != NULL) { reportSol(iter, change, TRUE, N, Z); }

    free(Zold);
  }

double pst_imgsys_sol_change(double* Zold, double* Z, int N)
  {
    int k;
    double max_change = 0;
    for(k = 0; k< N; k++)
      { double dk = Z[k] - Zold[k];
        if(dk < 0) { dk = - dk; }
        if(dk > max_change) { max_change = dk; }
      }
    return max_change;
  }

int* pst_imgsys_sort_equations(pst_imgsys_t *S)
  {
    int M = 2*(MAXCOEFS-1); /* Max number of arcs out of or into a node. */
    int N = S->N; /* Variables and equations are numbered from {0 to N-1}. */
    
    /* BUILDING THE PRIORITY GRAPH */

    /* The priority graph has one node per system variable,
      and at most one arc per non-diagonal coefficient of the system, 
      pointing from the unknown with higher priority to that
      of lower priority. It may have duplicate arcs. */
    
    /* Arcs of the priority graph: */
    int *dst = (int *)notnull(malloc(M*N*sizeof(int)), "no mem");
    /* Value of {n_out[k]} is count of neighbours of node {k} with lower priority. */
    int *n_out = (int *)notnull(malloc(N*sizeof(int)), "no mem");
    /* Value of {n_in[k]} is count of unprocessed neighbours of node {k} with higher priority. */
    int *n_in = (int *)notnull(malloc(N*sizeof(int)), "no mem");
    /* The arcs out of node {k} go to nodes {dst[M*k+j]} for {j} in {0..n_out[k]-1}. */

    auto void maybe_add_arrow(int k1, int k2);
      /* Adds an arc from node {k1} to node {k2}, 
        if the weight of {k1} is smaller than that of {k2}. */
    
    void maybe_add_arrow(int k1, int k2)
      { /* Nodes with lower weights have higher priority: */
        double w1 = S->eq[k1].wtot;
        double w2 = S->eq[k2].wtot;
        if (w1 < w2) { dst[M*k1 + n_out[k1]] = k2; n_out[k1]++;  n_in[k2]++; }
      }

    /* Gather the graph arcs: */
    int k;
    for(k = 0; k < N; k++){ n_in[k] = 0; n_out[k] = 0; }
    for(k = 0; k < N; k++)
      { pst_imgsys_equation_t *eqk = &(S->eq[k]);
        assert( eqk->ix[0] == k);
        int j;
        for(j = 1; j < eqk->nt; j++)
          { int i = eqk->ix[j];
            maybe_add_arrow(k,i);
            maybe_add_arrow(i,k);
          }
      }

    /* TOPOLOGICAL SORT OF THE PRIORITY GRAPH */
    
    /* A topological sort of the graph gives a list {ord[0..N-1]} 
      of nodes such that if {k1} has higher priority than {k2}
      then {k1} appears before {k2} in {ord}. */

    /* Queue of source nodes, and output ordering: */
    int *ord = (int *)notnull(malloc(N*sizeof(int)), "no mem"); /* Nodes in toporder. */
    int q_free = 0;  /* Index of first free entry in {ord}. */
    int q_start = 0; /* Index in {ord} of first unprocessed entry. */ 
    /* Nodes {ord[0..q_start-1]} have been sorted and have no outgoing or incoming edges. */
    /* Nodes {ord[q_start..q_free-1]} have been sorted and have no incoming edges. */

    auto void queue_insert(int index);
    void queue_insert(int index)
      { assert(q_free < N);
        ord[q_free] = index;
        q_free++;
      }

    auto int queue_remove(void);
    int queue_remove(void)
      { assert(q_start < q_free );
        int i = ord[q_start];
        q_start++;
        return i;
      }  

    /* Insert the nodes of highest priority: */
    for(k = 0; k < N; k++) { if (n_in[k] == 0) { queue_insert(k); } }

    while(q_start < q_free )
      { /* Get the next unprocessed node {k}: */
        int k = queue_remove();
        /* Add all its dominated neighbors to the queue, erase the arcs: */
        int j;
        for(j = 0; j < n_out[k]; j++)
          { /* Get a dominated neighbor {i}: */
            int i = dst[M*k + j];
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
    
void pst_imgsys_write(FILE *wr, pst_imgsys_t *S)
  {
    int k;
    bool_t fail = FALSE;
    for(k = 0; k < S->N; k++)
      {
        /* Get equation {k}: */
        pst_imgsys_equation_t *eqk = &(S->eq[k]);
        
        /* Compute indices of pixel corresponding to the unknown {Z[k]}: */
        int xk = S->col[k]; int yk = S->row[k];
         
        /* Print indices of main pixel: */
        fprintf(wr, "eq[%d][%d]:", xk, yk);

        /* Print the equation coefficients: */
        int t;
        for(t = 0; t < eqk->nt; t++)
          { /* Get hold of another unknown {Z[j]} that occurs in equation {k}: */
            int j = eqk->ix[t];
            /* Compute indices of the corresponding pixel. */
            int xj = S->col[j]; int yj = S->row[j];
            double cf_j = eqk->cf[t];
            if (t > 0) { fprintf(wr, " + "); }
            fprintf(wr, "%f*Z[%d][%d]", cf_j, xj, yj);
            if ((j < 0) || (j >= S->N))
              { fprintf(stderr, "pst_imgsys_write: invalid variable index %d in equation %d\n", j, k);
                fail = TRUE;
              }
          }
        /* Print the independent term: */
        fprintf(wr, " = %f", S->eq[k].rhs);
        fprintf(wr, "\n");
      }
    fflush(wr);
    demand(! fail, "write failed");
  }

void pst_imgsys_write_report(pst_imgsys_t *S, char *filePrefix, int level, char *tag, int indent)
  {
    if (S == NULL) { return; }
    char *fileName = float_image_mscale_file_name(filePrefix, level, -1, tag, "sys");
    fprintf(stderr, "%*swriting %s ...", indent, "", fileName);
    FILE* wr = fopen(fileName, "wt");
    assert(wr != NULL); 
    pst_imgsys_write(wr, S);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "\n");
    free(fileName);
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
