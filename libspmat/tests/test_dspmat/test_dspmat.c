#define PROG_NAME "test_dspmat"
#define PROG_DESC "test of {dspmat.h} (and therefore of {spmat.h})"
#define PROG_VERS "1.0"

#define test_dspmat_C_COPYRIGHT "Copyright © 2008  by the State University of Campinas (UNICAMP)"
/* Created on 2008-07-05 by J. Stolfi, UNICAMP */
/* Last edited on 2011-06-06 17:51:49 by stolfi */ 

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <jsrandom.h>
#include <jsfile.h>
#include <jsmath.h>
#include <affirm.h>
#include <bool.h>
#include <fget.h>
#include <affirm.h>


#include <spmat.h>
#include <spmat_io.h>
#include <spmat_linalg.h>
#include <dspmat.h>

#define i32min(X,Y) (((X) <= (Y) ? (X) : (Y)))
#define i32max(X,Y) (((X) >= (Y) ? (X) : (Y)))

/* PROTOTYPES */

void test_dspmat(int nt);
void test_dspmat_new(int it, bool_t verbose);  
void test_dspmat_add_element(int it, bool_t verbose);
void test_dspmat_trim(int it, bool_t verbose);
void test_dspmat_sort_entries(int it, bool_t verbose);
void test_dspmat_copy(int it, bool_t verbose);
void test_dspmat_sort_entries_ins(int it, bool_t verbose);
void test_dspmat_condense(int it, bool_t verbose);
void test_dspmat_mix(int it, bool_t verbose);
void test_dspmat_transpose(int it, bool_t verbose);
void test_dspmat_write_dspmat_read(int it, bool_t verbose);
void test_dspmat_extract_row_dspmat_extract_col_dspmat_mul(int it, bool_t verbose);
void test_dspmat_map_row(int it, bool_t verbose);
void test_dspmat_map_col(int it, bool_t verbose);
void test_dspmat_add_diagonal_dspmat_fill_diagonal(int it, bool_t fill, bool_t verbose);

void show_dspmat(char *Mname, dspmat_t *M, dspmat_count_t nPrint);

bool_t check_dspmat(dspmat_t *A, char *Aname, dspmat_t *R, char *Rname, bool_t verbose);
  /* Compares {A} with {R}; if {verbose}, also calls {show_dspmat} on them. */
  
void dspmat_throw(dspmat_t *M, double frac);
  /* Replaces the elements of {*M} by random elements.  Each element will be set,
    with  probability {frac}, to an independent random
    number distributed uniformly between 0 and 1;
    otherwise it will be set to the trivial value. */
  
void dspmat_scramble_entries(dspmat_t *M);
 /* Permutes the entries of {*M} in random order. */

void compare_vectors(double v[], char *vname, double r[], char *rname, dspmat_count_t n, double tol);
void check_order_of_entries(dspmat_t *A, int orow, int ocol);

/* IMPLEMENTATIONS */

int main (int argn, char **argv)
  { test_dspmat(30);  
    return 0;
  }

void test_dspmat(int nt)
  { fprintf(stderr, "Checking {dspmat_t} and its operations...\n");
    int it;
    for (it = 0; it < nt; it++)
      { 
        fprintf(stderr, "=== pass %d ===\n", it);
        bool_t verbose = (it < 4);
        test_dspmat_new(it, verbose);  
        test_dspmat_add_element(it, verbose);
        test_dspmat_trim(it, verbose);
        test_dspmat_sort_entries(it, verbose);
        test_dspmat_copy(it, verbose);
        test_dspmat_sort_entries_ins(it, verbose);
        test_dspmat_condense(it, verbose);
        test_dspmat_mix(it, verbose);
        test_dspmat_transpose(it, verbose);
        test_dspmat_extract_row_dspmat_extract_col_dspmat_mul(it, verbose);
        test_dspmat_map_row(it, verbose);
        test_dspmat_map_col(it, verbose);
        test_dspmat_add_diagonal_dspmat_fill_diagonal(it, FALSE, verbose);
        test_dspmat_add_diagonal_dspmat_fill_diagonal(it, TRUE, verbose);
        if (it < 3) 
          { test_dspmat_write_dspmat_read(it, verbose);
          }
         /* !!! Test the other operations !!! */
        fprintf(stderr, "\n");
      }
  }
   
void test_dspmat_new(int it, bool_t verbose)  
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_new} ...\n");
    
    dspmat_size_t Arows = 13, Acols = 11;
    dspmat_t A = dspmat_new(Arows,Acols,63); if (verbose) show_dspmat("A", &A, 1);
    assert(A.rows == Arows);
    assert(A.cols == Acols);
    assert(A.ents == 63);
    assert(A.e != NULL);
    
    dspmat_size_t Brows = 19, Bcols = 31;
    dspmat_t B = dspmat_new(Brows,Bcols, 0); if (verbose) show_dspmat("B", &B, 1);
    assert(B.rows == Brows);
    assert(B.cols == Bcols);
    assert(B.ents == 0);
    assert(B.e == NULL);
    
    dspmat_size_t Rrows = 49, Rcols = 79;
    dspmat_t R = dspmat_new(Rrows,Rcols, 0); if (verbose) show_dspmat("R", &R, 1);
    assert(R.rows == Rrows);
    assert(R.cols == Rcols);
    assert(R.ents == 0);
    assert(R.e == NULL);
    free(A.e);
    free(B.e);
    free(R.e);
  }

void test_dspmat_add_element(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_add_element} ...\n");
    
    dspmat_count_t ngen = 100;
    dspmat_size_t Arows = 13, Acols = 11;
    dspmat_t A = dspmat_new(Arows,Acols,12);
    dspmat_t R = dspmat_new(0,0,ngen);
    
    int k;
    dspmat_index_t rowMax = 0, colMax = 0;
    dspmat_pos_t posA = 0, posR = 0;
    for (k = 0; k < 100; k++)
      { int i = (k + (k*k/3)) / 19;
        int j = (k + (k*k/3)) % 19;
        if (i > rowMax) { rowMax = i; }
        if (j > colMax) { colMax = j; }
        double Aij = 1 + sin(i)*sin(3*j);
        posA = dspmat_add_element(&A, posA, i, j, Aij);
        if (Aij != 0.0) 
          { R.e[posR] = (dspmat_entry_t){ .row = i, .col = j, .val = Aij };
            posR++; 
          }
      }
    if (verbose) fprintf(stderr, "rowMax = %u  colMax = %u\n", rowMax, colMax);
    
    if (verbose) show_dspmat("A", &A, 1);
    assert(A.rows == i32max(Arows, rowMax+1));
    assert(A.cols == i32max(Acols, colMax+1));
    assert(A.ents >= posA);
    
    if (verbose) show_dspmat("R", &R, 1);
    R.rows = A.rows;
    R.cols = A.cols;
    assert(R.ents >= posR);
    assert(posA == posR);
    for (k = 0; k < posR; k++)
      { assert(A.e[k].row == R.e[k].row);
        assert(A.e[k].col == R.e[k].col);
        assert(A.e[k].val == R.e[k].val);
      }
    free(A.e);
    free(R.e);
  }
  
void test_dspmat_trim(int it, bool_t verbose)
  {
   if (verbose) fprintf(stderr, "\ntesting {dspmat_trim} ...\n");
    
    dspmat_count_t ngen = 100;
    dspmat_t A = dspmat_new(17,21,15);
    
    int k;
    dspmat_pos_t posA = 0;
    for (k = 0; k < ngen; k++)
      { int i = (k + (k*k/3)) / 19;
        int j = (k + (k*k/3)) % 19;
        double Aij = drandom();
        posA = dspmat_add_element(&A, posA, i, j, Aij);
      }
    assert(posA <= A.ents);
    dspmat_trim(&A, posA); assert(A.ents == posA); if (verbose) show_dspmat("A", &A, 1);
    posA = posA/2;
    dspmat_trim(&A, posA); assert(A.ents == posA); if (verbose) show_dspmat("A", &A, 1);
    
    free(A.e);
  }

void test_dspmat_sort_entries(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_sort_entries} ...\n");
    
    dspmat_t A = dspmat_new(11,17,0); 
    dspmat_throw(&A, 0.25);
    dspmat_scramble_entries(&A);
    if (verbose) show_dspmat("A", &A, 10);
    
    auto void do_test(int orow, int ocol);
    
    void do_test(int orow, int ocol)
      { if (verbose) fprintf(stderr, "ordering: orow = %+2d ocol = %+2d\n", orow, ocol);
        dspmat_scramble_entries(&A); 
        dspmat_sort_entries(&A, orow, ocol);
        if (verbose) show_dspmat("A", &A, 10); 
        check_order_of_entries(&A, orow, ocol);
      }
      
    do_test(+2, +1);
    do_test(+1, -2);
    do_test(-1, 00);
    free(A.e);
  }
  
void test_dspmat_copy(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_copy} ...\n");
    
    dspmat_t A = dspmat_new(11,17,0); 
    dspmat_throw(&A, 0.25); dspmat_scramble_entries(&A);
    
    dspmat_t R = dspmat_new(25,200,0);
    
    dspmat_copy(&A, &R);
    check_dspmat(&A, "A", &R, "R", verbose);
    
    free(A.e);
    free(R.e);
  }
  
void test_dspmat_sort_entries_ins(int it, bool_t verbose)
  {
    dspmat_t A = dspmat_new(11,17,0); 
    dspmat_throw(&A, 0.25); 
    dspmat_scramble_entries(&A);
    
    
    if (verbose) fprintf(stderr, "\ntesting {dspmat_sort_entries_ins} ...\n");
    dspmat_t R = dspmat_new(17,312,11);     
    dspmat_copy(&A, &R);
    
    dspmat_sort_entries_ins(&A, +1, +2, 0,A.ents); if (verbose) show_dspmat("A", &A, 10);
    dspmat_sort_entries(&R, +1, +2); if (verbose) show_dspmat("R", &R, 10);
    check_dspmat(&A, "A", &R, "R", verbose);
    
    dspmat_sort_entries_ins(&A, +2, +1, 0,A.ents); if (verbose) show_dspmat("A", &A, 10);
    dspmat_sort_entries(&R, +2, +1); if (verbose) show_dspmat("R", &R, 10);
    check_dspmat(&A, "A", &R, "R", verbose);
  }
  
void test_dspmat_condense(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_condense} ...\n");
    
    dspmat_t A = dspmat_new(0,0,0);
    dspmat_t R = dspmat_new(0,0,0);
    
    int k, r;
    int rMax = 5; /* Max repetitions of each index pair. */
    dspmat_pos_t posA = 0, posR = 0;
    for (r = 0; r < rMax; r++)
      { /* Generate some set of distinct pairs {i,j}: */
        for (k = 0; k < 300; k++)
          { int i = (k + (k*k/3)) / 19;
            int j = (k + (k*k/3)) % 19;
            int rCount = ((i + 3*j) % rMax) + 1; /* Number of repetitions for this {i,j}. */
            if (r < rCount)
              { double Aij = r;
                posA = dspmat_add_element(&A, posA, i, j, Aij);
              }
            if (r == 0)
              { double Rij = rCount*(rCount-1)/2;
                posR = dspmat_add_element(&R, posR, i, j, Rij);
              }
          }
      }
    dspmat_trim(&A, posA); dspmat_sort_entries(&A, +2, +1); if (verbose) show_dspmat("A", &A, 10);
    dspmat_trim(&R, posR); dspmat_sort_entries(&R, +2, +1); if (verbose) show_dspmat("R", &R, 10);

    auto double dbl_add(dspmat_index_t row, dspmat_index_t col, double v0, double v1);
    
    double dbl_add(dspmat_index_t row, dspmat_index_t col, double v0, double v1)
      { return v0+v1; }
    
    dspmat_condense(&A, dbl_add);
    check_dspmat(&A, "A", &R, "R", verbose);
    
    free(A.e);
    free(R.e);
  }

void test_dspmat_mix(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_mix} ...\n");
    
    dspmat_t A = dspmat_new(0,0,0);
    dspmat_t B = dspmat_new(0,0,0);
    dspmat_t R = dspmat_new(0,0,0);
    
    int k;
    dspmat_pos_t posA = 0, posB = 0, posR = 0;
    for (k = 0; k < 100; k++)
      { int i = (k + (k*k/3)) / 19;
        int j = (k + (k*k/3)) % 19;
        double Aij = 1 + sin(i)*sin(3*j);
        double Bij = 1 + cos(i)*cos(5*j);
        double Rij = 0.5 * Aij + 2.0 * Bij;
        posA = dspmat_add_element(&A, posA, i, j, Aij);
        posB = dspmat_add_element(&B, posB, i, j, Bij);
        posR = dspmat_add_element(&R, posR, i, j, Rij);
      }
    dspmat_trim(&A, posA); dspmat_sort_entries(&A, +2, +1); if (verbose) show_dspmat("A", &A, 10);
    dspmat_trim(&B, posB); dspmat_sort_entries(&B, +2, +1); if (verbose) show_dspmat("B", &B, 10);
    dspmat_trim(&R, posR); dspmat_sort_entries(&R, +2, +1); if (verbose) show_dspmat("R", &R, 10);
    
    dspmat_t C = dspmat_new(0,0,0);
    dspmat_mix(0.5, &A, 2.0, &B, &C);
    check_dspmat(&C, "C", &R, "R", verbose);
    
    free(A.e);
    free(B.e);
    free(C.e);
    free(R.e);
  }
  
void test_dspmat_transpose(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_transpose} ...\n");
    
    dspmat_t A = dspmat_new(5,7,0); 
    dspmat_throw(&A, 0.25); 
    dspmat_scramble_entries(&A);
    
    dspmat_transpose(&A);  if (verbose) show_dspmat("A", &A, 10);
    free(A.e);
  }
  
void test_dspmat_write_dspmat_read(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_write} ...\n");
    
    dspmat_t A = dspmat_new(5,7,0); dspmat_throw(&A, 0.25); dspmat_scramble_entries(&A);
    FILE *wr = open_write("out/dspmat.txt", TRUE);
    dspmat_write(wr, &A);
    fclose(wr);
    
    if (verbose) fprintf(stderr, "\ntesting {dspmat_read} ...\n");
    dspmat_t R = dspmat_new(0,0,0);
    FILE *rd = open_read("out/dspmat.txt", TRUE);
    dspmat_read(rd, &R);
    fclose(rd);
    
    check_dspmat(&A, "A", &R, "R", verbose);
    free(A.e);
    free(R.e);
  }

void test_dspmat_extract_row_dspmat_extract_col_dspmat_mul(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat,dspmat_extract_col,dspmat_mul} ...\n");
    
    dspmat_t A = dspmat_new(7,11,0); 
    dspmat_throw(&A, 0.50); 
    if (verbose) show_dspmat("A", &A, 10);
    
    dspmat_t B = dspmat_new(11,5,0); 
    dspmat_throw(&B, 0.50); 
    if (verbose) show_dspmat("B", &B, 10);
    
    dspmat_t C = dspmat_new(0,0,0);
    /* Multiply by the library, result is {C}: */
    dspmat_mul(&A, &B, &C); if (verbose) show_dspmat("C", &C, 10);
    
    /* Multiply by extracting rows and cols, result is {R}: */
    double Arow[A.cols];
    double Bcol[B.rows];
    dspmat_sort_entries(&A, +2, +1) /* (by row,col) */; if (verbose) show_dspmat("A", &A, 10);
    dspmat_sort_entries(&B, +1, +2) /* (by col,row) */; if (verbose) show_dspmat("B", &B, 10);
    
    dspmat_t R = dspmat_new(0,0,0);
    dspmat_pos_t posR = 0;
    int i;
    dspmat_pos_t posA = 0;
    for (i = 0; i < A.rows; i++)
      { posA = dspmat_extract_row(&A, posA, i, Arow, A.cols);
        int j;
        dspmat_pos_t posB = 0; 
        for (j = 0; j < B.cols; j++)
          { posB = dspmat_extract_col(&B, posB, j, Bcol, B.rows);
            double sum = 0;
            int k;
            for (k = 0; k < A.cols; k++) { sum += Arow[k]*Bcol[k]; }
            posR = dspmat_add_element(&R, posR, i,j, sum);
          }
      }
    dspmat_trim(&R, posR);
    
    check_dspmat(&C, "C", &R, "R", verbose);
    free(A.e);
    free(B.e);
    free(C.e);
    free(R.e);
  }
  
void test_dspmat_map_row(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_map_row} ...\n");
    
    dspmat_t A = dspmat_new(7,11,0);
    dspmat_throw(&A, 0.50);
    if (verbose) show_dspmat("A", &A, 10);
    
    dspmat_size_t nu = A.rows;
    double u[nu];
    dspmat_size_t nv = A.cols;
    double v[nv];
    int i,j;
    for (i = 0; i < nu; i++) { u[i] = drandom(); }
    /* By {dspmat_map_row}, result in {v}: */
    dspmat_map_row(u, nu, &A, v, nv);
    /* By the book, reult in {r}: */
    dspmat_size_t nr = A.cols;
    double r[nr];
    dspmat_sort_entries_ins(&A, +1, +2, 0, A.ents) /* (by col,row) */; if (verbose) show_dspmat("A", &A, 10);
    double Acol[A.rows];
    dspmat_pos_t posA = 0;
    for (j = 0; j < A.cols; j++) 
      { double sum = 0; 
        posA = dspmat_extract_col(&A, posA, j, Acol, A.rows);
        assert(nu == A.rows);
        for (i = 0; i < nu; i++) { sum += u[i]*Acol[i]; }
        r[j]= sum;
      }
    assert(posA == A.ents);
    compare_vectors(v, "v", r, "r", nv, 1.0e-10);
    free(A.e);
  }
  
void test_dspmat_map_col(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_map_col} ...\n");
    
    dspmat_t A = dspmat_new(7,11,0);
    dspmat_throw(&A, 0.50);
    if (verbose) show_dspmat("A", &A, 10);
    
    dspmat_size_t nu = A.cols;
    double u[nu];
    dspmat_size_t nv = A.rows;
    double v[nv];
    int i,j;
    for (i = 0; i < nu; i++) { u[i] = drandom(); }
    /* By {dspmat_map_col}, result in {v}: */
    dspmat_map_col(&A, u, nu, v, nv);
    /* By the book, reult in {r}: */
    dspmat_size_t nr = A.rows;
    double r[nr];
    dspmat_sort_entries(&A, +2, +1) /* (by row, col) */; if (verbose) show_dspmat("A", &A, 10);
    double Arow[A.cols];
    dspmat_pos_t posA = 0;
    for (i = 0; i < A.rows; i++) 
      { double sum = 0; 
        posA = dspmat_extract_row(&A, posA, i, Arow, A.cols);
        assert(A.cols = nu);
        for (j = 0; j < nu; j++) { sum += Arow[j]*u[j]; }
        r[i]= sum;
      }
    compare_vectors(v, "v", r, "r", nv, 1.0e-10);
    free(A.e);    
  }
  
void test_dspmat_add_diagonal_dspmat_fill_diagonal(int it, bool_t fill, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "\ntesting {dspmat_add_diagonal,dspmat_fill_diagonal} ...\n");
    
    dspmat_t C = dspmat_new(17,31,0);
    dspmat_count_t nd = imin(lcm(C.cols, C.rows), 300);
    double d[nd];
    int k;
    double vfill = drandom();
    for (k = 0; k < nd; k++) { d[k] = (fill ? vfill : drandom()); }
    dspmat_pos_t posC = 0;
    dspmat_index_t row1 = 5, col1 = 3;
    if (fill)
      { posC = dspmat_fill_diagonal(&C, posC, row1, col1, vfill, nd); }
    else
      { posC = dspmat_add_diagonal(&C, posC, row1, col1, d, nd); }
    dspmat_trim(&C, posC);
    
    /* Check elements: */
    dspmat_t R = dspmat_new(C.rows,C.cols,0);
    dspmat_pos_t posR = 0;
    for (k = 0; k < nd; k++) 
      { dspmat_index_t row = (row1 + k) % R.rows;
        dspmat_index_t col = (col1 + k) % R.cols;
        posR = dspmat_add_element(&R, posR, row, col, d[k]);
      }
    dspmat_trim(&R, posR);
         
    check_dspmat(&C, "C", &R, "R", verbose);
    free(C.e);   
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
      
void compare_vectors(double v[], char *vname, double r[], char *rname, dspmat_size_t n, double tol)
  { int k;
    for (k = 0; k < n; k++) 
      { if ((v[k] != r[k]) && (fabs(rel_diff(v[k], r[k])) > tol))
          { fprintf(stderr, "** vector element mismatch\n");
            fprintf(stderr, "  %s[%d] = %24.16e\n", vname, k, v[k]);
            fprintf(stderr, "  %s[%d] = %24.16e\n", rname, k, r[k]);
            fprintf(stderr, "  error = %24.16e\n", v[k] - r[k]);
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

void dspmat_scramble_entries(dspmat_t *M)
  {
    dspmat_pos_t p;
    for (p = 0; p < M->ents; p++)
      { dspmat_pos_t q = p + (int)floor(drandom()*(M->ents - p));
        if (q >= M->ents) { q = M->ents - 1; }
        dspmat_entry_t e = M->e[p];
        M->e[p] = M->e[q];
        M->e[q] = e;
      }
  }
  
void check_order_of_entries(dspmat_t *A, int orow, int ocol)
{
  dspmat_pos_t p1;
  for (p1 = 1; p1 < A->ents; p1++)
    { dspmat_pos_t p0 = p1 - 1;
      dspmat_entry_t *e0 = &(A->e[p0]);
      dspmat_entry_t *e1 = &(A->e[p1]);
      assert((e0->col < A->cols) && (e0->row < A->rows));
      assert((e1->col < A->cols) && (e1->row < A->rows));
      int cmp = spmat_compare_indices(e0->row, e0->col, e1->row, e1->col, orow, ocol);
      if (cmp > 0) 
        { fprintf(stderr, "** entries out of order (args = %+d %+d)\n", orow, ocol);
          fprintf(stderr, "  A->e[%d] = [%d][%d] (%24.16e)\n", p0, e0->row, e0->col, e0->val);
          fprintf(stderr, "  A->e[%d] = [%d][%d] (%24.16e)\n", p1, e1->row, e1->col, e1->val);
          assert(FALSE);
        }
    }
}



