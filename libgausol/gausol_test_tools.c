/* see gausol_test_tools.h  */
/* Last edited on 2024-11-30 04:53:40 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include <jsrandom.h>
#include <bool.h>
#include <affirm.h>

#include <gausol_print.h>
#include <gausol_test_tools.h>

/* INTERNAL PROTOTYPES */

double gausol_test_tools_max_det_roundoff(uint32_t m, double rms_A);
  /* Estimates the maximum roundoff error in the computation of the
    determinant of an {m × m} matrix, before or after the triangulation,
    given the RMS value {rms_A} of the elements of the original matrix
    and the cleanup threshold {tiny} given in the triangulation. */

/* The do { .. } while is a hak to allow semicolon after it. */
#define BAD(fmt, ...) \
  do { ok = FALSE; fprintf(stderr, "** " fmt "\n", ##__VA_ARGS__); } while (FALSE)

/* IMPLEMENTATIONS */

void gausol_test_tools_choose_system
  ( uint32_t trial,
    uint32_t m_max, uint32_t n_max, uint32_t p_max,
    uint32_t *m_P, 
    uint32_t *n_P, double **A_P,
    uint32_t *p_P, double **B_P,
    double **X_P,
    double tiny,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }
    srand(1665 + 12*trial);
    srandom(1665 + 12*trial);

    uint32_t sptrials = 4;  /* Num trial with special systems */
    uint32_t sptrials_mult = 4; /* How many time to repeat the same special trial. */

    uint32_t m, n, p;
    if (trial < sptrials*sptrials_mult)
      { m = trial/sptrials_mult+1; 
        n = trial/sptrials_mult+1; 
        p = 2;
      }
    else
      { m = uint32_abrandom(1, m_max); /* Rows (equations). */
        n = uint32_abrandom(1, n_max); /* Main columns (unknowns). */
        p = uint32_abrandom(1, p_max); /* RHS columns (problems). */
      }

    if (verbose) { fprintf(stderr, "    choosing m = %d  n = %d  p = %d\n", m, n, p); }

    double *A = NULL; /* Main systems matrix. */
    double *X; /* Nominal solution. */

    if (verbose) { fprintf(stderr, "    generating system...\n\n"); }
    if (trial < sptrials*sptrials_mult)
      { /* A fixed system for each trial: */
        A = talloc(m*n, double);
        X = talloc(n*p, double);
        for (uint32_t j = 0; j < n; j++)
          { for (uint32_t i = 0; i < m; i++)
              { double Aij = n*sin((i+1)*(j+1) + trial + M_SQRT2);
                A[i*n + j] = Aij;
              }
            for (uint32_t k = 0; k < p; k++)
              { double Xjk = m*cos((j+1)*(k+1) + trial + M_LN2);
                X[j*p + k] = Xjk;
              }
          }
      }
    else
      { /* A random system: */
        gausol_test_tools_throw_system(m, n, &A, p, &X, verbose);
      }

    /* Compute the right-hand side {B = A X}: */
    double *B = talloc(m*p, double);  /* Righ-hand side. */
    gausol_test_tools_multiply(m, n, p, A, X, B, tiny);
    if (verbose) 
      { gausol_print_system
          ( stderr, 6, "%12.6f", "original system:", 
            m,NULL,0,
            n,NULL,0, "A",A, 
            p,"B",B,  
            0,NULL,NULL, 
            ""
          );
        gausol_print_array
          ( stderr, 6, "%12.6f", "nominal solution", 
            n,NULL,0, p,NULL,0, "X", X, ""
          );
      }

    (*m_P) = m; (*n_P) = n; (*p_P) = p;
    (*A_P) = A; (*B_P) = B; (*X_P) = X;
    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
  }
    
void gausol_test_tools_multiply
  ( uint32_t m, uint32_t n, uint32_t p,
    double A[], double X[], double B[],
    double tiny
  )
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    if (debug) { fprintf(stderr, "      m = %d  n = %d  p = %d  tiny = %24.16e\n", m, n, p, tiny); }
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t k = 0; k < p; k++)
          { double sum = 0.0, corr = 0.0;
            for (uint32_t j = 0; j < n; j++)
              { double Aij = A[i*n + j];
                double Xjk = X[j*p + k];
                double term = Aij*Xjk;
                /* Kahan's summation: */
                double tcorr = term - corr;
                double newSum = sum + tcorr;
                corr = (newSum - sum) - tcorr;
                sum = newSum;
              }
            if (fabs(sum) < tiny) 
              { if (debug) { fprintf(stderr, "    sum = %24.16e cleared\n", sum); }
                sum = 0.0;
              }
            B[i*p + k] = sum;
          }
      }
      if (debug) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
  }

void gausol_test_tools_check_triang_reduce
  ( uint32_t m, uint32_t prow[], 
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    double tiny,
    uint32_t rank, 
    double det_ref, double det_cmp,
    double rms_A,
    bool_t verbose
  )
  { 
    if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }
    /* Check the {gausol_triang_reduce} post-conditions: */
    demand((rank <= m) && (rank <= n), "invalid rank {rank}");
    if (prow != NULL) { for (uint32_t i = 0; i < m; i++) { demand(prow[i] < m, "invalid index in {prow}"); } }
    if (pcol != NULL) { for (uint32_t j = 0; j < n; j++) { demand(pcol[j] < n, "invalid index in {pcol}"); } }

    if (verbose) 
      { gausol_print_system
          ( stderr, 6, "%12.6f", "triangulated system", 
            m,prow,0,  n,pcol,0, "A",A,  p,"B",B,  0,NULL,NULL, ""
          ); 
      }
    
    /* Check that {P00} is upper diagonal and {P10} is zero: */
    bool_t ok = TRUE;
    for (uint32_t t = 0; t < rank; t++)
      { uint32_t prow_t = (prow == NULL ? t : prow[t]); assert(prow_t < m);
        uint32_t pcol_t = (pcol == NULL ? t : pcol[t]); assert(pcol_t < n);
        double Ptt = A[prow_t*n + pcol_t];
        if ((! isfinite(Ptt)) || (fabs(Ptt) < 1.0e-180))
          { BAD("diagonal element P[%d,%d] = %24.16e is invalid or too small", t, t, Ptt); }
        for (uint32_t i = t+1; i < m; i++)
          { uint32_t prow_i = (prow == NULL ? i : prow[i]); assert(prow_i < m);
            double Pit = A[prow_i*n + pcol_t];
            if ((! isfinite(Pit)) || (Pit != 0)) 
              { BAD("sub-diagonal element P[%d,%d] = %24.16e is invalid or not zero", i, t, Pit); }
          }
      }
    
    /* Check that the relevant part of {P11} satisfies (3.**): */
    uint32_t ilim = ((prow == NULL) && (rank < m) ? rank + 1 : m);
    uint32_t jlim = ((pcol == NULL) && (rank < n) ? rank + 1 : n); 
    for (uint32_t i = rank; i < ilim; i++)
      { uint32_t prow_i = (prow == NULL ? i : prow[i]); assert(prow_i < m);
        for (uint32_t j = rank; j < jlim; j++)
          { uint32_t pcol_j = (pcol == NULL ? j : pcol[j]); assert(pcol_j < n);
            double Pij = A[prow_i*n + pcol_j];
            if (Pij != 0) { BAD("supposedly zero element P[%d,%d] = %24.16e is not zero", i, j, Pij); }
          }
      } 
    demand(ok, "{gausol_triang_reduce} failed");
    
    /* Check determinant: */
    if (! isnan(det_ref))
      { if (verbose) { fprintf(stderr, "  checking preservation of determinant ...\n"); }
        assert(m == n);
        gausol_test_tools_compare_determinants( m, n, rank, rms_A, det_ref, det_cmp, verbose);
      }

    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
  }

void gausol_test_tools_check_diagonalize
  ( uint32_t m, uint32_t prow[], 
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    uint32_t rank,
    bool_t verbose
  )
  { if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }

    if (verbose) 
      { gausol_print_system
          ( stderr, 6, "%12.6f", "diagonalized system", 
            m,prow,0,  n,pcol,0, "A",A,  p,"B",B,  0,NULL,NULL, ""
          ); 
      }
    
    bool_t ok = TRUE;
    for (uint32_t i = 0; i < rank; i++)
      { uint32_t prow_i = (prow == NULL ? i : prow[i]); assert(prow_i < m);
        for (uint32_t j = i + 1; j < rank; j++)
          { uint32_t pcol_j = (pcol == NULL ? j : pcol[j]); assert(pcol_j < n);
            double *Pij = &(A[prow_i*n + pcol_j]);
            if ((*Pij) != 0) { BAD("nonzero element P[%d,%d] = %24.16e above diagonal\n", i, j, (*Pij)); }
          }
      }
    demand(ok, "** matrix is not diagonal");
    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
  }
  
void gausol_test_tools_check_normalize
  ( uint32_t m, uint32_t prow[],
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    uint32_t rank,
    bool_t verbose
  )
  { if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }

    if (verbose) 
      { gausol_print_system
          ( stderr, 6, "%12.6f", "normalized system", 
            m,prow,0,  n,pcol,0, "A",A,  p,"B",B,  0,NULL,NULL, ""
          ); 
      }

    bool_t ok = TRUE;
    for (uint32_t t = 0; t < rank; t++)
      { uint32_t prow_t = (prow == NULL ? t : prow[t]); assert(prow_t < m);
        uint32_t pcol_t = (pcol == NULL ? t : pcol[t]); assert(pcol_t < n);
        double Ptt = A[prow_t*n + pcol_t];
        if (Ptt != 1.0) { BAD("diagonal element A[%d,%d] = %24.16e is not 1\n", t, t, Ptt); }
      }
    demand(ok, "** matrix is not normalized");
    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
  }

void gausol_test_tools_compare_determinants
  ( uint32_t m, uint32_t n,
    uint32_t rank,
    double rms_A,
    double det_ref, double det_cmp,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }
    demand(m == n, "main matrix should be square");
    demand(rank <= m, "invalid rank {rank}");
    if (! (isnan(det_ref) || isnan(det_cmp)))
      { double tol = gausol_test_tools_max_det_roundoff(m, rms_A);
        if (verbose) { fprintf(stderr, "    determinant: enum =  %24.16e triang = %24.16e\n", det_ref, det_cmp); }
        demand(isfinite(det_ref), "overflow in reference determinant");
        demand(isfinite(det_cmp), "overflow in determinant from triangulation");
        double dif = det_cmp - det_ref;
        if (verbose || (fabs(dif) > tol)) 
          { fprintf(stderr, "    determinant difference =  %24.16e tol = %24.16e\n", dif, tol); 
            demand(fabs(dif) <= tol, "** determinants do not match");
          }
      }
    else 
      { if (isnan(det_ref))
          { if (verbose) { fprintf(stderr, "    the reference determinant was not computed"); } }
        if (isnan(det_cmp))
          { if (verbose) { fprintf(stderr, "    {gausol_triang_reduce} could not compute the determinant"); } }
      }
    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
  }

void gausol_test_tools_check_satisfaction
  ( uint32_t m,  
    uint32_t n,  double A[],
    uint32_t p, double B[],
    double X_ref[],
    bool_t verbose
  )
  { if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }
    double rel_tol = 1.0e-14;
    double abs_tol = 1.0e-15;
    gausol_test_tools_check_residual
      ( m, n, A, p, B, X_ref, rel_tol, abs_tol, verbose);
    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
  }

double gausol_test_tools_check_solve
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p,  double B[],
    double X_ref[], double X_cmp[],
    uint32_t rank,
    bool_t verbose
  )
  { 
    if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }

    if (verbose) 
      { gausol_print_array
          ( stderr, 6, "%12.6f", "computed solution", 
            n,NULL,0,  p,NULL,0, "X_cmp",X_cmp,  ""
          ); 
      }
      
    double rms_X_ref = gausol_test_tools_elem_RMS(n, p, X_ref);
    double rms_X_cmp = gausol_test_tools_elem_RMS(n, p, X_cmp);
    double rms_X = hypot(rms_X_ref, rms_X_cmp)/M_SQRT2;

    double rel_tol = 3.0e-14;
    double abs_tol = 3.0e-14*rms_X;
    double err_X = NAN;
    demand((rank <= m) && (rank <= n), "invalid rank {rank}");
    if ((rank == m) && (rank == n))
      { /* The solution should be unique, so: */
        err_X = gausol_test_tools_compare_solutions
          (n, p, X_ref, X_cmp, rel_tol, abs_tol, verbose);
      }
    if (rank == n)
      { /* The problem is solvable, so: */
        gausol_test_tools_check_residual
          (m, n, A, p, B, X_cmp, rel_tol, abs_tol, verbose);
      }
    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
    return err_X;
  }  

double gausol_test_tools_compare_solutions
  ( uint32_t n, uint32_t p,
    double X_ref[], double X_cmp[],
    double rel_tol, double abs_tol,
    bool_t verbose
  )
  { if (verbose) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    bool_t ok = TRUE;
    double sum_e2 = 0;
    for (uint32_t k = 0;  k < p; k++)
      { for (uint32_t j = 0;  j < n; j++)
          { double Xjk_ref = X_ref[j*p + k];
            double Xjk_cmp = X_cmp[j*p + k];
            demand(! isnan(Xjk_cmp), "invalid {X_cmp} element");
            demand(fabs(Xjk_cmp) < 1.0e+180, "{X_cmp} element too large");
            double dif = Xjk_cmp - Xjk_ref;
            sum_e2 += dif*dif;
            double mag = (fabs(Xjk_ref) + fabs(Xjk_cmp))/2;
            double err_max = hypot(abs_tol, rel_tol*mag);
            if (fabs(dif) > err_max)
              { BAD("solution mismatch (X_ref-X_cmp)[%d,%d] = %24.16e max = %24.16e\n", j, k, dif, err_max); }
          }
      }
    demand(ok, "** computed solution does not match reference solution");
    if (verbose) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
    return sqrt(sum_e2/(n*p));
  }

void gausol_test_tools_check_residual
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p, double B[],
    double X[],
    double rel_tol, double abs_tol,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    /* Computes the max abs elems in each row of {A} and each col of {X}: */
    double Amax[m];
    double Xmax[p];
    for (uint32_t i = 0; i < m; i++) { Amax[i] = 0.0; }
    for (uint32_t k = 0; k < p; k++) { Xmax[k] = 0.0; }
    for (uint32_t j = 0; j < n; j++)
      { for (uint32_t i = 0; i < m; i++) 
          { Amax[i] = fmax(Amax[i], fabs(A[i*n + j])); }
        for (uint32_t k = 0; k < p; k++) 
          { Xmax[k] = fmax(Xmax[k], fabs(X[j*p + k])); }
      }
      
    /* Compute {B_cmp = A X}: */
    double *B_cmp = talloc(m*p, double);
    double Btiny = 0.0;
    gausol_test_tools_multiply(m, n, p, A, X, B_cmp, Btiny);
      
    /* Compare {B_cmp} with the original {B}: */
    bool_t ok = TRUE;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t k = 0;  k < p; k++)
          { double Bik_ref = B[i*p + k];
            double Bik_cmp = B_cmp[i*p + k];
            double err = Bik_cmp - Bik_ref;
            double err_max = hypot(abs_tol, rel_tol*Amax[i]*Xmax[k]);
            if (fabs(err) > err_max)
              { BAD("residual mismatch (A X - B)[%d,%d] = %24.16e max = %24.16e\n", i, k, err, err_max); }
          }
      }
    demand(ok, "** residual {A X - B} is not small enough");
    free(B_cmp);
    if (verbose) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
  }

void gausol_test_tools_make_row_dependent(uint32_t i, uint32_t m, uint32_t n, double A[])
  { 
    demand(i < m, "invalid row index {i}");
    for (uint32_t j = 0; j < n; j++) { A[i*n + j] = 0.0; }
    for (uint32_t k = 0; k < m; k++)
      { if (k != i)
          { /* Add a random multiple of row {k} to row {i}: */
            double alpha = dabrandom(-1.0, +1.0);
            for (uint32_t j = 0; j < n; j++) 
              { A[i*n + j] += alpha*A[k*n + j]; }
          }
      }
  }
  
double gausol_test_tools_det_by_enum(uint32_t m, uint32_t n, double A[], uint32_t q)
  { /* !!! Duplicates {rmxn_det_by_enum} from {libgeo}. Cleanup somehow. !!! */
    if ((q > m) || (q > n)) { return 0.0; }
    demand(q <= gausol_test_tools_det_by_enum_SIZE_MAX, "too many determinant terms to enumerate");
    
    uint32_t perm[q];

    for (uint32_t i = 0;  i < q; i++) { perm[i] = i; }

    double det = 0.0, corr = 0.0;

    auto void add_terms(uint32_t k, double product);
      /* Generates all permutations of {perm[k..q-1]}, and adds to
        {det} the partial {product} times {A[rank,perm[rank]} for {rank} in {k..q-1}.
        Upon return, {perm[k..q-1]} are back to their initial
        state. */

    add_terms(0, 1.0);
    return det;

    void add_terms(uint32_t k, double product)
      { if (k >= q)
          { /* The {product} is a complete determinant term */
            double term = product; 
            /* Kahan's detmation: */
            double tcorr = term - corr;
            double newDet = det + tcorr;
            corr = (newDet - det) - tcorr;
            det = newDet;
          }
        else
          { for (uint32_t j = k; j < q; j++)
              { if (j != k) { uint32_t t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
                uint32_t row = k, col = perm[k];
                double prod_more = product * (k == j ? 1.0 : -1.0) * A[row*n + col];
                add_terms(k+1, prod_more);
                if (j != k) { uint32_t t = perm[k]; perm[k] = perm[j]; perm[j] = t; }
              }
          }
      }
  }

double gausol_test_tools_elem_RMS(uint32_t m, uint32_t n, double A[])
  { double sum = 0, corr = 0;
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { double el = A[i*n + j];
            double term = el*el; 
            /* Kahan's detmation: */
            double tcorr = term - corr;
            double newSum = sum + tcorr;
            corr = (newSum - sum) - tcorr;
            sum = newSum;
          }
      }
    return sqrt(sum/(m*n));
  }

double *gausol_test_tools_throw_matrix
  ( uint32_t m, uint32_t n,
    double pzero,
    double mag_min, double mag_max,
    double tiny,
    char *head, 
    char *name,
    bool_t verbose
  )
  {
    demand((pzero >= 0) && (pzero <= 1.0), "invalid {pzero}");
    demand(mag_min > 1.0e-15, "invalid {mag_min}");
    demand((mag_max >= mag_min) && (mag_max < 1.0e+15), "invalid {max_max}");
    double mag_range = mag_max/mag_min;
    
    double *A = talloc(m*n, double);
    if ((m > 0) && (n > 0))
      { double pelz = (m < n ? pzero/m : pzero/n); /* Probab of elem being zero. */
        for (uint32_t i = 0;  i < m; i++)
          { for (uint32_t j = 0; j < n; j++)
             { double el;
               if (drandom() < pelz)
                 { el = 0.0; }
               else
                 { /* Generate {mag} uniform in log scale: */
                   double mag = mag_min * pow(mag_range, drandom());
                   /* Generate {el} uniform in {[-mag _ +mag]}: */
                   el = dabrandom(-mag, +mag);
                   if (fabs(el) < tiny) { el = 0; }
                 }
               A[i*n + j] = el;
             }
          }
      }
    if (verbose)
      { gausol_print_array
         ( stderr, 6, "%12.6f", head, 
           m,NULL,0,  n,NULL,0, name,A, ""
         );
      }
    return A;
  }

void gausol_test_tools_throw_system
  ( uint32_t m,
    uint32_t n, double **A_P,
    uint32_t p, double **X_P,
    bool_t verbose
  )
  {
    /* Generate power-of-ten scale factors: */
    double mag_A = pow(10.0, int32_abrandom(-3, +3));
    double mag_X = pow(10.0, int32_abrandom(-3, +3));
    if (verbose) { fprintf(stderr, "  scales: A = %8.1e X = %8.1e\n", mag_A, mag_X); }
    
    double tiny = 1.0e-15; /* Min abs value of nonzero elem. */
    double pzero = 0.25;   /* Avg count of zeros in each row or col. */

    /* Generate a random coefficient matrix {A}: */
    (*A_P) = gausol_test_tools_throw_matrix(m, n, pzero, mag_A, mag_A, tiny, "original coefficient matrix", "A", verbose);

    /* Generate a random solution matrix {X_ref}: */
    (*X_P) = gausol_test_tools_throw_matrix(n, p, pzero, mag_X, mag_X, tiny, "original solution", "X", verbose);
  }

double gausol_test_tools_max_det_roundoff(uint32_t m, double rms_A)
  { double nf = m; /* Number of factors in each term */
    double nt = 1.0; /* Number of terms. */
    for (uint32_t i = 1; i <= m; i++) { nt *= i; }
    double rel_tiny = 1.0e-14;
    /* Assumes that the roundoff errors are independent: */
    return rel_tiny * pow(rms_A, nf) * sqrt(nt);
  }

