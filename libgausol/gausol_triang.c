/* See gausol_triang_reduce.h */
/* Last edited on 2024-11-30 15:03:59 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>

#include <gausol_triang.h>

void gausol_triang_reduce
  ( uint32_t m, uint32_t prow[],
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    double tiny,
    double *det_P,
    uint32_t *rank_P
  )
  { 
    bool_t debug = FALSE;
    double sgn = +1;
    uint32_t rank = 0;
    
    auto bool_t pivot(void);
      /* Assumes {rank<m} and {rank<n}. Looks for the best pivot in the
        relevant part of submatrix {P11}.
        Returns {FALSE} if and only the pivot candidates are all zero.
        Otherwise updates {pcol[rank..n-1]} and/or
        {prow[rank..m-1]} so that the pivot becomes {P11[0,0] = P[rank,rank]}.
        Inverts the sign of {sgn} if these changes required an odd
        number of pair swaps. 
        
        If oth {prow} and {pcol} are {NULL}, considers only {P11[0,0]},
        that is, {P[rank,rank]}, as the candidate for pivot.
        
        If {prow} is not {NULL} but {pcol} is {NULL}, looks only in the
        first column of {P11}.   May change {prow[rank..m-1]}. 
        
        If {prow} is {NULL} but {pcol} is not {NULL}, looks only in the
        first row of {P11}. May change {pcol[rank..n-1]}.
        
        If both {prow} and{pcol} are non-{NULL}, looks for the pivot in
        the whole sumbatrix {P11}. May change both {prow[rank..m-1]} and
        {pcol[rank..n-1]}.
        
        In any case, if there are two or more candidates, chooses the
        one that has maximum absolute value.
        
        In any case, if this procedure returns {FALSE}, the searched
        part of {P11} is all zero, meaning that one or more of
        conditions (3.00), (3.10), (3.01), or (3.11) are satisfied,
        according to whether {prow} and/or {pcol} are {NULL} or not. */
       
    auto void clear_eqs(void);
      /* Assumes that {rank<m} and {rank<n}, conditions (1) and (2) are satisfied,
        and {P[rank,rank]} is nonzero.
        
        For each row {i} in {rank+1..m-1}, computes {alpha = P[i,rank]/P[rank,rank]},
        and sets {P[i,rank]} to zero. Then, for each column {j} in {rank+1..n-1},
        sets {P[i,j] = P[i,j] - alpha*P[rank,j]}. For each {k} in {0..p-1},
        sets {Q[i,k] = Q[i,j] - alpha*Q[rank,k]}.  If any of those
        modified elements become less than {tiny} in absolute value,
        sets that element to zero. */
        
    auto double determinant(void);
      /* Returns the product of the diagonal elements of the permuted
        submatrix {P00}, times {sgn}. If the matrix is square, {rank} is
        {m}, this will be the determinant of {A}. */

    if ((m == 0) && (n == 0)) 
      { /* Nothing to do. */ }
    else
      { /* Initialize {pcol} and/or {prow} as appropriate: */
        if ((m != 0) && (prow != NULL)) { for (uint32_t i = 0; i < m; i++) { prow[i] = i; } }
        if ((n != 0) && (pcol != NULL)) { for (uint32_t j = 0; j < n; j++) { pcol[j] = j; } }

        /* Main loop: */
        while ((rank < m) && (rank < n))
          { bool_t found_pivot = pivot();
            if (! found_pivot) { break; }
            clear_eqs();
            rank++;
          }
      }
      
    if (det_P != NULL)
      { (*det_P) = determinant(); }
   
    if (rank_P != NULL) { (*rank_P) = rank; }
    return;
    
    bool_t pivot(void)
      { if (debug) { fprintf(stderr, "    enter {pivot} rank = %d\n", rank); }
        assert((rank < m) && (rank < n));
        double qual_best = 0;  /* Best pivot quality indicator so far. */
        uint32_t i_best, j_best; /* Row and col in {P} of best pivot so far, if {qual_best > 0}. */
        uint32_t ilim = ((prow == NULL) && (rank < m) ? rank+1 : m);
        uint32_t jlim = ((pcol == NULL) && (rank < n) ? rank+1 : n);
        for (uint32_t i = rank; i < ilim; i++)
          { uint32_t prow_i = (prow == NULL ? i : prow[i]);
            for (uint32_t j = rank; j < jlim; j++)
              { uint32_t pcol_j = (pcol == NULL ? j : pcol[j]); 
                double *Aij = &(A[prow_i*n + pcol_j]);
                if (rank == 0)
                  { /* First look at element, check if OK: */ 
                    demand(isfinite(*Aij) && (fabs(*Aij) <= 1.0e+180), "invalid matrix element");
                    if (fabs(*Aij) < tiny) { (*Aij) = 0; } 
                  }
                double qual = fabs(*Aij);
                if (qual > qual_best)
                  { i_best = i; j_best = j; qual_best = qual; }
              }
          }
        if (qual_best == 0)
          { if (debug) { fprintf(stderr, "      no pivot found\n"); }
            return FALSE; 
          }
        else
          { if (debug) { fprintf(stderr, "      found pivot at P[%d,%d] qual = %24.16e\n", i_best, j_best, qual_best); }
            /* Adjust {prow} if needed: */
            if (i_best != rank) 
              { assert(prow != NULL);
                uint32_t t = prow[i_best]; prow[i_best] = prow[rank]; prow[rank] = t;
                sgn = -sgn;
              }
            /* Adjust {pcol} if needed: */
            if (j_best != rank)
              { assert(pcol != NULL);
                uint32_t t = pcol[j_best]; pcol[j_best] = pcol[rank]; pcol[rank] = t;
                sgn = -sgn;
              }
            return TRUE;
          }
      }
      
    void clear_eqs(void)
      { if (debug) { fprintf(stderr, "    enter {clear_eqs} rank = %d\n", rank); }
        assert((rank < m) && (rank < n));
        uint32_t prow_r = (prow == NULL ? rank : prow[rank]); assert(prow_r < m);
        uint32_t pcol_r = (pcol == NULL ? rank : pcol[rank]); assert(pcol_r < n);
        double *Prr = &(A[prow_r*n + pcol_r]); /* The pivot element. */
        if (debug) { fprintf(stderr, "      pivot P[%d,%d]=A[%d,%d] = %24.16e\n", rank, rank, prow_r, pcol_r, (*Prr)); }
        assert(fabs(*Prr) >= tiny);
        for (uint32_t i = rank+1; i < m; i++)
          { uint32_t prow_i = (prow == NULL ? i : prow[i]); assert(prow_i < m);
            double *Pir = &(A[prow_i*n + pcol_r]);
            double alpha = (*Pir)/(*Prr);
            if (debug) { fprintf(stderr, "      clearing row P[%d,*],Q[%d,*] = A[%d,*],B[%d,*] alpha = %24.16e\n", i, i, prow_i, prow_i, alpha); }
            demand(isfinite(alpha), "overflow in row operation (1)"); 
            (*Pir) = 0.0;
            for (uint32_t j = rank+1; j < n; j++)
              { uint32_t pcol_j = (pcol == NULL ? j : pcol[j]); assert(pcol_j < n);
                double *Pij = &(A[prow_i*n + pcol_j]);
                double *Prj = &(A[prow_r*n + pcol_j]);
                double Pij_new = (*Pij) - alpha*(*Prj);
                demand(isfinite(Pij_new), "overflow in row operation (2)"); 
                if (fabs(Pij_new) < tiny) { Pij_new = 0.0; }
                (*Pij) = Pij_new;
              }
            for (uint32_t k = 0; k < p; k++)
              { double *Qik = &(B[prow_i*p + k]);
                double *Qrk = &(B[prow_r*p + k]);
                double Qik_new = (*Qik) - alpha*(*Qrk);
                demand(isfinite(Qik_new), "overflow in row operation (3)"); 
                if (fabs(Qik_new) < tiny) { Qik_new = 0.0; }
                (*Qik) = Qik_new;
              }
           }
       }
       
    double determinant(void)
      { double det = sgn;
        for (uint32_t t = 0; t < rank; t++)
          { uint32_t prow_t = (prow == NULL ? t : prow[t]); assert(prow_t < m);
            uint32_t pcol_t = (pcol == NULL ? t : pcol[t]); assert(pcol_t < n);
            double Ptt = A[prow_t*n + pcol_t];
            demand(isfinite(Ptt), "triangulation produced invalid diag elem");
            det *= Ptt;
          }
        demand(isfinite(det), "overflow in determinant computation");
        return det;
      }
  }

void gausol_triang_diagonalize
  ( uint32_t m, uint32_t prow[],
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    uint32_t rank,
    double tiny
  )
  { 
    demand((rank <= m) && (rank <= n), "invalid rank {rank}");
    for (uint32_t t = 1; t < rank; t++)
      { uint32_t prow_t = (prow == NULL ? t : prow[t]); assert(prow_t < m);
        uint32_t pcol_t = (pcol == NULL ? t : pcol[t]); assert(pcol_t < n);
        double Ptt = A[prow_t*n + pcol_t];
        demand(isfinite(Ptt) && (fabs(Ptt) >= 1.0e-180), "invalid {P00} diagonal element");
        for (uint32_t i = 0; i < t; i++)
          { uint32_t prow_i = (prow == NULL ? i : prow[i]); assert(prow_i < m);
            double *Pit = &(A[prow_i*n + pcol_t]);
            double alpha = (*Pit)/Ptt;
            if (fabs(*Pit) >= tiny)
              { (*Pit) = 0;
                for (uint32_t j = t+1; j < n; j++)
                  { uint32_t pcol_j = ((pcol == NULL) ? j : pcol[j]);  assert(pcol_j < n);
                    double *Pij = &(A[prow_i*n + pcol_j]);
                    if (fabs(*Pij) < tiny) 
                      { demand((*Pij) == 0, "tiny elem of {P00,P01} was not cleared (1)"); }
                    double *Ptj = &(A[prow_t*n + pcol_j]);
                    if (fabs(*Ptj) < tiny) 
                      { demand((*Ptj) == 0, "tiny elem of {P00,P01} was not cleared (2)"); }
                    (*Pij) -= alpha*(*Ptj);
                    if (fabs(*Pij) < tiny) { (*Pij) = 0; }
                  }
                for (uint32_t k = 0; k < p; k++)
                  { double *Qik = &(B[prow_i*p + k]);
                    if (fabs(*Qik) < tiny) 
                      { demand((*Qik) == 0, "tiny elem of {Q0} was not cleared (1)"); }
                    double *Qtk = &(B[prow_t*p + k]);
                    if (fabs(*Qtk) < tiny) 
                      { demand((*Qtk) == 0, "tiny elem of {Q0} was not cleared (2)"); }
                    (*Qik) -= alpha*(*Qtk);
                    if (fabs(*Qik) < tiny) { (*Qik) = 0; }
                  }
                
              }
            else
              { demand((*Pit) == 0, "tiny element of {P00} was not cleared"); }
          }
      }
  }

void gausol_triang_normalize
  ( uint32_t m, uint32_t prow[],
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    uint32_t rank,
    double tiny
  )
  { demand(rank <= m, "invalid rank {rank}");
    for (uint32_t t = 0; t < rank; t++)
      { uint32_t prow_t = (prow == NULL ? t : prow[t]); assert(prow_t < m);
        uint32_t pcol_t = (pcol == NULL ? t : pcol[t]); assert(pcol_t < n);
        double *Ptt = &(A[prow_t*n + pcol_t]);
        demand(isfinite(*Ptt) && (fabs(*Ptt) >= 1.0e-180), "invalid {P00} diagonal element");
        double piv = (*Ptt);
        (*Ptt) = 1.0;
        for (uint32_t j = rank; j < n; j++)
          { uint32_t pcol_j = ((pcol == NULL) ? j : pcol[j]); assert(pcol_j < n);
            double *Ptj = &(A[prow_t*n + pcol_j]);
            (*Ptj) /= piv;
            if (fabs(*Ptj) < tiny) { (*Ptj) = 0; }
          }
        for (uint32_t k = 0; k < p; k++)
          { double *Qtk = &(B[prow_t*p + k]);
            (*Qtk) /= piv;
            if (fabs(*Qtk) < tiny) { (*Qtk) = 0; }
          }
      }
  }

