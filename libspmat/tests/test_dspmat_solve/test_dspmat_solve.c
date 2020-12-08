#define PROG_NAME "test_dspmat_solve"
#define PROG_DESC "test of lineat algebra modules for {dspmat.h}"
#define PROG_VERS "1.0"

#define test_dspmat_solve_C_COPYRIGHT "Copyright © 2007  by the State University of Campinas (UNICAMP)"
/* Created on 2007-01-02 by J. Stolfi, UNICAMP */
/* Last edited on 2018-03-04 23:00:22 by stolfilocal */ 

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <bool.h>
#include <affirm.h>

#include <dspmat.h>
#include <dspmat_extra.h>
#include <dspmat_linsys_GS.h>
#include <dspmat_linsys_ALT.h>
/* !!! #include <dspmat_linsys_BOOT.h> !!! */

typedef enum 
  { solver_GS,
    solver_ALT,
    solver_BOOT
  } solver_t;

/* PROTOTYPES */

#define i32min(X,Y) (((X) <= (Y) ? (X) : (Y)))
#define i32max(X,Y) (((X) >= (Y) ? (X) : (Y)))

void test_dspmat_extra(int nt);
void test_dspmat_solve(int it, bool_t verbose, solver_t solver);

void show_dspmat(char *Mname, dspmat_t *M, dspmat_count_t nPrint);
  /* Prints the first {nPrint} entries of {M} to {stderr}. */
  
void show_vec(char *vname, double v[], dspmat_size_t n, dspmat_size_t nPrint);
  /* Prints the first {nPrint} entries of {v[0..n-1]} to {stderr}. */

bool_t check_dspmat(dspmat_t *A, char *Aname, dspmat_t *R, char *Rname, bool_t verbose);
  /* Compares {A} with {R}; if {verbose}, also calls {show_dspmat} on them. */
  
void dspmat_throw(dspmat_t *M, double frac);
  /* Replaces the elements of {*M} by random elements.  Each element will be set,
    with  probability {frac}, to an independent random
    number distributed uniformly between 0 and 1;
    otherwise it will be set to the trivial value. */
  
void dspmat_throw_nzd(dspmat_t *M, double frac, double mag);
  /* Same as {dspmat_throw} but makes sure that the matrix has
    diagonal elements with absolute value {mag} or greater. */
  
void compare_vectors(double v[], char *vname, double r[], char *rname, dspmat_size_t n, double abs_tol, double rel_tol);

/* IMPLEMENTATIONS */

int main (int argn, char **argv)
  { test_dspmat_extra(30);  
    return 0;
  }

void test_dspmat_extra(int nt)
  { fprintf(stderr, "Checking {dspmat_extra,dspmat_linsys_{GS,ALT,BOOT}} ...\n");
    int it;
    for (it = 0; it < nt; it++)
      { 
        fprintf(stderr, "=== pass %d ===\n", it);
        bool_t verbose = (it < 4);
        test_dspmat_solve(it, verbose, solver_GS);  
        test_dspmat_solve(it, verbose, solver_ALT);
        fprintf(stderr, "** skipping the test of the BOOT solver\n");
        /* !!! test_dspmat_solve(it, verbose, solver_BOOT); !!! */
        /* !!! Test the other operations !!! */
        fprintf(stderr, "\n");
      }
  }
   
void test_dspmat_solve(int it, bool_t verbose, solver_t solver)
  {
    char *solver_name[3] = {"GS", "ALT", "BOOT"};
    if (verbose) 
      { fprintf(stderr, "\ntesting {dspmat_linsys_%s} ...\n", solver_name[solver]); }
    
    /* Generate a random linear equation system {A,b} of random size {n}: */
    dspmat_size_t n = int32_abrandom(1,11)*int32_abrandom(1,17);
    dspmat_t A = dspmat_new(n,n,0); 
    double b[n];
    int i;
    if (it == 0)
      { /* Use a diagonal matrix and a simple right-hand-side vector: */ 
        dspmat_pos_t posA = dspmat_fill_diagonal(&A, 0, 0,0, 4.0, n);
        dspmat_trim(&A, posA);
        for (i = 0; i < A.rows; i++) { b[i] = 8+i; }
      }
    else
      { /* Use a random matrix with nonzero diagonal, and a random vector: */  
        dspmat_throw_nzd(&A, 0.20, 1.0); 
        for (i = 0; i < n; i++) { b[i] = 2*(drandom() - 0.5); }
      }
    if (verbose) show_dspmat("A", &A, 10);
    if (verbose) show_vec("b", b, n, 10);
    
    /* Solve it: */
    double x[n];
    int max_iter = 30;
    double omega = 0.50;
    double abs_tol = 1.0e-4;
    double rel_tol = 1.0e-3;
    switch (solver)
      { 
        case solver_GS:
          dspmat_linsys_GS_solve(b, n, &A, x, n, max_iter, omega, abs_tol, rel_tol);
          break;
        case solver_ALT: 
          dspmat_linsys_ALT_solve(b, n, &A, x, n, max_iter, omega, abs_tol, rel_tol);
          break;
        case solver_BOOT: 
          /* !!! dspmat_linsys_BOOT_solve(b, n, &A, x, n, max_iter, omega, abs_tol, rel_tol); !!! */
          break;
        default:
          assert(FALSE);
      }
    
    /* Check solution: */
    double y[n];  /* Result of {A*x}. */
    dspmat_map_col(&A, x, n, y, n);
    compare_vectors(b, "b", y, "y", n, abs_tol, rel_tol);
    
    dspmat_trim(&A,0);
  }

bool_t check_dspmat(dspmat_t *A, char *Aname, dspmat_t *R, char *Rname, bool_t verbose)
  {
    if (verbose)
      { fprintf(stderr, "check_dspmat\n");
        show_dspmat(Aname, A, 1);
        show_dspmat(Rname, R, 1);
        /* !!! To be completed !!! */
      }
    if (((A->rows) != (R->rows)) || ((A->cols) != (R->cols)) || ((A->ents) != (R->ents)))
      { fprintf(stderr, "** size mismatch:\n");
        fprintf(stderr, "  %s = %5d rows  %5d cols  %5d ents\n", Aname, A->rows, A->cols, A->ents);
        fprintf(stderr, "  %s = %5d rows  %5d cols  %5d ents\n", Rname, R->rows, R->cols, R->ents);
        assert(FALSE);
      }
    
    dspmat_pos_t p;
    for (p = 0; p < A->ents; p++)
      { dspmat_entry_t *a = &(A->e[p]);
        dspmat_entry_t *r = &(R->e[p]);
        if ((a->row != r->row) || (a->col != r->col) || (a->val != r->val))
          { fprintf(stderr, "** entries don't match\n");
            fprintf(stderr, "  %s.e[%d] = [%d][%d] (%24.16e)\n", Aname, p, a->row, a->col, a->val);
            fprintf(stderr, "  %s.e[%d] = [%d][%d] (%24.16e)\n", Rname, p, r->row, r->col, r->val);
            assert(FALSE);
          }
      }
    return TRUE;
  }
      
void compare_vectors(double v[], char *vname, double r[], char *rname, dspmat_size_t n, double abs_tol, double rel_tol)
  { int k;
    for (k = 0; k < n; k++) 
      { double vk = v[k];
        double rk = r[k];
        double dbase2 = abs_tol*abs_tol + rel_tol*rel_tol*vk*vk + 1.0e-300;
        double dk = vk - rk;
        double error = dk/sqrt(dbase2);
        if ((v[k] != r[k]) && (fabs(error) > 1.0))
          { fprintf(stderr, "** vector element mismatch\n");
            fprintf(stderr, "  %s[%d] = %24.16e\n", vname, k, v[k]);
            fprintf(stderr, "  %s[%d] = %24.16e\n", rname, k, r[k]);
            fprintf(stderr, "  diff   = %24.16e\n", dk);
            fprintf(stderr, "  error  = %17.15e\n", error);
            assert(FALSE);
          }
      }
  }

void show_dspmat(char *Mname, dspmat_t *M, dspmat_count_t nPrint)
  {
    if (nPrint > M->ents) { nPrint = M->ents; }
    fprintf(stderr, "%s = { %3u %3u %10p[%5u] } = (", Mname, M->rows, M->cols, M->e, M->ents); 
    if (M->ents > 0)
      { dspmat_pos_t k;
        for (k = 0; k < nPrint; k++)
          { dspmat_entry_t *eP = &(M->e[k]); 
            if (nPrint > 1) { fprintf(stderr, "\n   "); }
            fprintf(stderr, " [%3u][%3u] = %+24.16e", eP->row, eP->col, eP->val);
          }
        
        if (M->ents > nPrint) 
          { if (nPrint > 1) { fprintf(stderr, "\n   "); }
            fprintf(stderr, " ...");
          }
        
      }
    if (nPrint > 1) { fprintf(stderr, "\n "); }
    fprintf(stderr, " )\n");
  }
   
void dspmat_throw(dspmat_t *M, double frac)
  { dspmat_index_t i,j;
    dspmat_pos_t posM = 0;
    for (i = 0; i < M->rows; i++)
      { for (j = 0; j < M->cols; j++)
          { double Mij = (drandom() < frac ? drandom() : 0.0);
            posM = dspmat_add_element(M, posM, i, j, Mij);
          }
      }
    dspmat_trim(M, posM);
  }

void dspmat_throw_nzd(dspmat_t *M, double frac, double mag)
  { 
    dspmat_throw(M, frac);
    
    /* Make sure that entries are sorted: */
    dspmat_sort_entries(M, +2, +1);
    
    /* Make sure that diagonal elements are at leas {mag}: */
    dspmat_count_t n_old = M->ents;
    dspmat_pos_t pos_old = 0;;
    dspmat_pos_t pos_new = M->ents; /* Current number of filled entries in {M}. */
    
    int row;
    for (row = 0; row < M->rows; row++)
      { bool_t has_diag = FALSE; /* TRUE iff row {row} has a diagonal elem. */
        while (pos_old < n_old)
          { dspmat_entry_t *mP = &(M->e[pos_old]);
            assert(mP->row >= row);
            if (mP->row > row) { break; }
            if (mP->col == row)
              { if (fabs(mP->val) < mag) { mP->val = copysign(mag, mP->val); } 
                has_diag = TRUE;
              }
            pos_old++;
          }
        if (! has_diag)
          { pos_new = dspmat_add_element(M, pos_new, row, row, copysign(mag, drandom()-0.5)); }
      }
    dspmat_trim(M, pos_new);
    
    /* Make sure that entries are sorted: */
    dspmat_sort_entries(M, +2, +1);
  }

void show_vec(char *vname, double v[], dspmat_size_t n, dspmat_size_t nPrint)
  { int i;
    for (i = 0; i < n; i++)
      { fprintf(stderr, "  %s[%4d] = %9.4f  %24.16e\n", vname, i, v[i], v[i]); }
  }


