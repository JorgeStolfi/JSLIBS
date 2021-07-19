#define PROG_NAME "test_float_array"
#define PROG_DESC "test of {float_array.h} (and therefore of {array.h})"
#define PROG_VERS "1.0"

#define test_float_array_C_COPYRIGHT "Copyright © 2009  by the State University of Campinas (UNICAMP)"
/* Created on 2009-08-31 by J. Stolfi, UNICAMP */
/* Last edited on 2021-07-18 19:16:49 by jstolfi */ 

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


#include <ix.h>
#include <array.h>
#include <array_io.h>
// #include <array_linalg.h>
#include <float_array.h>

#define i32min(X,Y) (((X) <= (Y) ? (X) : (Y)))
#define i32max(X,Y) (((X) >= (Y) ? (X) : (Y)))

#define szFMT ix_index_t_FMT

/* PROTOTYPES */

void test_float_array(int nt); 

float_array_t *float_array_new_random(int it, bool_t verbose);
  /* Creates an array with random number of axes, random size vector, 
    and random element order, using {float_array_new}. Also does 
    some consistency checks on the result. */

void test_float_array_new(int it, bool_t verbose);
void test_float_array_get_elem_set_elem(int it, bool_t verbose);
void test_float_array_copy( int it, bool_t verbose);
void test_float_array_new_descr(int it, bool_t verbose);
void test_float_array_copy_descr(int it, bool_t verbose);
void test_float_array_free_elems(int it, bool_t verbose);
void test_float_array_free_descr(int it, bool_t verbose);
void test_float_array_get_size(int it, bool_t verbose);
void test_float_array_check_size(int it, bool_t verbose);
void test_float_array_write_float_array_read(int it, bool_t verbose);

void show_float_array(FILE *wr, char *Mname, float_array_t *A, float_array_count_t nPrint);
void show_indices(FILE *wr, char *pf, ix_dim_t d, ix_index_t ix[], char *sf);

bool_t check_float_array(float_array_t *A, char *Aname, float_array_t *R, char *Rname, bool_t verbose);
  /* Compares {A} with {R}; if {verbose}, also calls {show_float_array} on them. */
  
void float_array_throw(float_array_t *A, double frac);
  /* Replaces the elements of {*A} by random elements.  Each element will be set,
    with  probability {frac}, to an independent random
    number distributed uniformly between 0 and 1;
    otherwise it will be set to the trivial value. */
  
void float_array_scramble_entries(float_array_t *A);
 /* Permutes the entries of {*A} in random order. */

void compare_vectors(double v[], char *vname, double r[], char *rname, float_array_count_t n, double tol);
void check_order_of_entries(float_array_t *A, int orow, int ocol);

/* IMPLEMENTATIONS */

int main (int argn, char **argv)
  { test_float_array(30);  
    return 0;
  }

void test_float_array(int nt)
  { fprintf(stderr, "Checking {float_array_t} and its operations...\n");
    int it;
    for (it = 0; it < nt; it++)
      { 
        fprintf(stderr, "=== pass %d ===\n", it);
        bool_t verbose = (it < 4);
        test_float_array_new(it, verbose);
        test_float_array_get_elem_set_elem(it, verbose);
        // test_float_array_copy(it, verbose);
        // test_float_array_copy_descr(it, verbose);
        // test_float_array_free_elems(it, verbose);
        // test_float_array_free_descr(it, verbose);
        // test_float_array_get_elem_position(it, verbose);
        // test_float_array_get_size(it, verbose);
        // test_float_array_check_size(it, verbose);
        if (it < 7) 
          { test_float_array_write_float_array_read(it, verbose); }
         /* !!! Test the other operations !!! */
        fprintf(stderr, "\n");
      }
  }
   
float_array_t *float_array_new_random(int it, bool_t verbose)
  {
    float_array_dim_t na = (float_array_dim_t)(it % array_MAX_AXES);
    float_array_size_t sz[na];
    int ia;
    ix_count_t ne = 1;
    for (ia = 0; ia < na; ia++) 
      { sz[ia] = int32_abrandom(0, 3)+int32_abrandom(0, 4); ne *= (ix_count_t)(sz[ia]); }
    ix_order_t ixor = (ix_order_t)int32_abrandom(0, 1);
    if (verbose)
      { fprintf(stderr, "%s(%d, %c)\n", __FUNCTION__, it, "FT"[verbose]);
        fprintf(stderr, "na = %d\n", na);
        fprintf(stderr, "sz = [");
        for (ia = 0; ia < na; ia++) { fprintf(stderr, " " szFMT, sz[ia]); }
        fprintf(stderr, " ]\n");
        fprintf(stderr, "ixor = %c\n", "LF"[ixor]);
      }
    float_array_t *A = float_array_new_descr();
    (*A) = float_array_new(na, sz);
    ix_descr_t *DA = &(A->ds);
    if (verbose) { show_float_array(stderr, "A", A, 5); }
    assert(DA->na == na);
    for (ia = 0; ia < na; ia++) { assert(DA->sz[ia] == sz[ia]); }
    if (ne == 0)
      { assert(A->e == NULL); }
    else
      { assert(A->e != NULL); }
    return A;
  }

void test_float_array_new(int it, bool_t verbose)  
  {
    if (verbose) fprintf(stderr, "Checking {float_array_new} ...\n");
    srandom(4634 + 17*(2*it + 1));
    
    float_array_t *A = float_array_new_random (it, verbose);
    
    if (! ix_descr_is_valid(&(A->ds), FALSE)) 
      { fprintf(stderr, "** invalid descriptor\n");
        show_float_array(stderr, "A", A, 5);
        assert(FALSE);
      }
    float_array_free_elems(A); float_array_free_descr(A);
  }

void test_float_array_get_elem_set_elem(int it, bool_t verbose)
  { if (verbose) fprintf(stderr, "Checking {float_array_{get_elem,get_elem_pos,set_elem}} ...\n");
    srandom(4634 + 19*(2*it + 1));
    
    float_array_t *A = float_array_new_random (it, verbose);
    ix_descr_t *DA = &(A->ds);
    float_array_dim_t na = DA->na;
    ix_index_t ix[na];
    
    /* If the array has replication along any axes, make them trivial, */
    /* otherwise the test becomes much more difficult. */
    int ia;
    for (ia = 0; ia < na; ia++)
      { if ((DA->sz[ia] > 1) && (DA->st[ia] == 0)) { DA->sz[ia] = 1; } }
    
    /* Fill all elements with distinct numbers: */
    int k = 4615;
    if (ix_descr_indices_first(DA, ix))
      { do {
          float v_set = (float)k;
          float_array_set_elem(A, ix, v_set);
          k++;
        } while (! ix_descr_next(DA,ix,NULL));
      }
    
    /* Check whether the numbers are still there: */
    k = 4615;
    if (ix_descr_indices_first(DA, ix))
      { ix_pos_t p_exp = ix_descr_position(DA, ix);
        do {
          float v_set = (float)k;
          float v_get = float_array_get_elem(A, ix);
          if (v_get != v_set)
            { fprintf(stderr, "** set/get elem mismatch:\n");
              show_float_array(stderr, "A", A, 5);
              show_indices(stderr, "ix = ", DA->na, ix, "\n");
              fprintf(stderr, "  v (set) = %15.8e  v (get) = %15.8e\n", v_set, v_get);
              assert(FALSE);
            }
          ix_pos_t p_get = float_array_get_elem_position(A, ix);
          if ((p_get != p_exp) || (A->e[p_get] != v_set))
            { fprintf(stderr, "** elem position mismatch:\n");
              show_float_array(stderr, "A", A, 5);
              show_indices(stderr, "ix = ", DA->na, ix, "\n");
              fprintf(stderr, "  p (get) = %12lu  p (exp) = %12lu\n", p_get, p_exp);
              fprintf(stderr, "  v (set) = %15.8e  v (get) = %15.8e\n", v_set, A->e[p_get]);
              assert(FALSE);
            }
          k++;
        } while (! ix_descr_next(DA,ix,&p_exp));
      }
    
  }

void test_float_array_write_float_array_read(int it, bool_t verbose)
  {
    if (verbose) fprintf(stderr, "Checking {float_array_write} ...\n");
    char *fname = NULL;
    asprintf(&fname, "out/float_array-%04d.txt", it);

    srandom(4634 + 23*(2*it + 1));
    float_array_t *A = float_array_new_random(it, verbose);
    
    float_array_throw(A, 0.25);
    FILE *wr = open_write(fname, TRUE);
    float_array_write(wr, A);
    fclose(wr);
    
    if (verbose) fprintf(stderr, "Checking {float_array_read} ...\n");
    FILE *rd = open_read(fname, TRUE);
    float_array_t *R = float_array_read(rd);
    fclose(rd);
    
    check_float_array(A, "A", R, "R", verbose);
    float_array_free_elems(A); float_array_free_descr(A);
    float_array_free_elems(R); float_array_free_descr(R);
    free(fname);
  }

bool_t check_float_array(float_array_t *A, char *Aname, float_array_t *R, char *Rname, bool_t verbose)
  {
    ix_descr_t *DA = &(A->ds);
    ix_descr_t *DR = &(R->ds);
    
    if (verbose)
      { fprintf(stderr, "check_float_array\n");
        show_float_array(stderr, Aname, A, 1);
        show_float_array(stderr, Rname, R, 1);
        /* !!! To be completed !!! */
      }
    if ((DA->na) != (DR->na))
      { fprintf(stderr, "** num axes mismatch:\n");
        fprintf(stderr, "  %s.ds.na = %5d ds.na\n", Aname, DA->na);
        fprintf(stderr, "  %s.ds.na = %5d ds.na\n", Rname, DR->na);
        assert(FALSE);
      }
      
    float_array_dim_t na = DA->na;
    int ia;
    for (ia = 0; ia < na; ia++)
      {
        if ((DA->sz[ia]) != (DR->sz[ia]))
          { fprintf(stderr, "** size mismatch in axis %d:\n", ia);
            fprintf(stderr, "  %s.ds.sz[%d] = %5ld\n", Aname, ia, DA->sz[ia]);
            fprintf(stderr, "  %s.ds.sz[%d] = %5ld\n", Rname, ia, DR->sz[ia]);
            assert(FALSE);
          }
      }
      
    if ((A->e == NULL) != (R->e == NULL))
      { fprintf(stderr, "** element pointer mismatch:\n");
        fprintf(stderr, "  %s.e = %10p\n", Aname, A->e);
        fprintf(stderr, "  %s.e = %10p\n", Rname, R->e);
        assert(FALSE);
      }
    
    ix_index_t ix[na];
    if (ix_descr_indices_first(DA, ix))
      { do {
          ix_pos_t pA = ix_descr_position(DA, ix);
          ix_pos_t pR = ix_descr_position(DR, ix);
          if (A->e[pA] != R->e[pR])
            { fprintf(stderr, "** element value mismatch\n");
              fprintf(stderr, "ix = [");
              for (ia = 0; ia < na; ia++) { fprintf(stderr, " %2ld", ix[ia]); }
              fprintf(stderr, " ]\n");
              fprintf(stderr, "  %s.e[ix] = %+15.8e\n", Aname, A->e[pA]);
              fprintf(stderr, "  %s.e[ix] = %+15.8e\n", Rname, R->e[pR]);
              assert(FALSE);
            }
        } while (! ix_descr_next(DA,ix,NULL));
      }
    return TRUE;
  }
      
void show_float_array(FILE *wr, char *Aname, float_array_t *A, float_array_count_t nPrint)
  {
    ix_descr_t *DA = &(A->ds);
    float_array_dim_t na = DA->na;
    ix_count_t ne = ix_descr_num_positions(DA);
    if (nPrint > ne) { nPrint = ne; }
    int ia;
    fprintf(wr, "%s = { na: %3u", Aname, na); 
    fprintf(wr, " sz: ["); 
    for (ia = 0; ia < na; ia++) { fprintf(wr, " %2ld", DA->sz[ia]); }
    fprintf(wr, " ]"); 
    fprintf(wr, " st: ["); 
    for (ia = 0; ia < na; ia++) { fprintf(wr, " %ld", DA->st[ia]); }
    fprintf(wr, " ]"); 
    fprintf(wr, " e: %10p[%5lu]", A->e, ne); 
    fprintf(wr, " } = ("); 
    float_array_count_t k = 0;
    ix_index_t ix[na];
    if (ix_descr_indices_first(DA, ix))
      { ix_pos_t pA = ix_descr_position(DA, ix);
        do {
          if (k >= nPrint)
            { fprintf(wr, " ...");
              break;
            }
          float v = A->e[pA];
          fprintf(wr, " %+15.8e", v);
        } while (! ix_descr_next(DA,ix,&pA));
      }
    fprintf(wr, " )\n");
  }
   

void show_indices(FILE *wr, char *pf, ix_dim_t d, ix_index_t ix[], char *sf)
  { fprintf(wr, "%s", pf); 
    int i;
    for (i = 0; i < d; i++)
      { fprintf(wr, " %2ld", ix[i]); }
    fprintf(wr, "%s", sf);
  }

void float_array_throw(float_array_t *A, double frac)
  { ix_descr_t *DA = &(A->ds);
    float_array_dim_t na = DA->na;
    ix_index_t ix[na];
    if (ix_descr_indices_first(DA, ix))
      { ix_pos_t pA = ix_descr_position(DA, ix);
        do {
          double v = (drandom() < frac ? drandom() : 0.0);
          A->e[pA] = (float)v;
        } while (! ix_descr_next(DA,ix,&pA));
      }
  }

// void compare_vectors(double v[], char *vname, double r[], char *rname, float_array_dim_t n, double tol)
//   { int k;
//     for (k = 0; k < n; k++) 
//       { if ((v[k] != r[k]) && (fabs(rel_diff(v[k], r[k])) > tol))
//           { fprintf(stderr, "** vector element mismatch\n");
//             fprintf(stderr, "  %s[%d] = %24.16e\n", vname, k, v[k]);
//             fprintf(stderr, "  %s[%d] = %24.16e\n", rname, k, r[k]);
//             fprintf(stderr, "  error = %24.16e\n", v[k] - r[k]);
//             assert(FALSE);
//           }
//       }
//   }
// 
// void test_float_array_add_elem(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_add_elem} ...\n");
//     
//     float_array_count_t ngen = 100;
//     float_array_dim_t Arows = 13, Acols = 11;
//     float_array_t A = float_array_new(Arows,Acols,12);
//     float_array_t R = float_array_new(0,0,ngen);
//     
//     int k;
//     float_array_index_t rowMax = 0, colMax = 0;
//     float_array_pos_t posA = 0, posR = 0;
//     for (k = 0; k < 100; k++)
//       { int i = (k + (k*k/3)) / 19;
//         int j = (k + (k*k/3)) % 19;
//         if (i > rowMax) { rowMax = i; }
//         if (j > colMax) { colMax = j; }
//         double Aij = 1 + sin(i)*sin(3*j);
//         posA = float_array_add_elem(&A, posA, i, j, Aij);
//         if (Aij != 0.0) 
//           { R.e[posR] = (float_array_entry_t){ .row = i, .col = j, .val = Aij };
//             posR++; 
//           }
//       }
//     if (verbose) fprintf(stderr, "rowMax = %u  colMax = %u\n", rowMax, colMax);
//     
//     if (verbose) show_float_array(stderr, "A", &A, 1);
//     assert(A.rows == i32max(Arows, rowMax+1));
//     assert(A.cols == i32max(Acols, colMax+1));
//     assert(A.ents >= posA);
//     
//     if (verbose) show_float_array(stderr, "R", &R, 1);
//     R.rows = A.rows;
//     R.cols = A.cols;
//     assert(R.ents >= posR);
//     assert(posA == posR);
//     for (k = 0; k < posR; k++)
//       { assert(A.e[k].row == R.e[k].row);
//         assert(A.e[k].col == R.e[k].col);
//         assert(A.e[k].val == R.e[k].val);
//       }
//     free(A.e);
//     free(R.e);
//   }
//   
// void test_float_array_trim(int it, bool_t verbose)
//   {
//    if (verbose) fprintf(stderr, "Checking {float_array_trim} ...\n");
//     
//     float_array_count_t ngen = 100;
//     float_array_t A = float_array_new(17,21,15);
//     
//     int k;
//     float_array_pos_t posA = 0;
//     for (k = 0; k < ngen; k++)
//       { int i = (k + (k*k/3)) / 19;
//         int j = (k + (k*k/3)) % 19;
//         double Aij = drandom();
//         posA = float_array_add_elem(&A, posA, i, j, Aij);
//       }
//     assert(posA <= A.ents);
//     float_array_trim(&A, posA); assert(A.ents == posA); if (verbose) show_float_array(stderr, "A", &A, 1);
//     posA = posA/2;
//     float_array_trim(&A, posA); assert(A.ents == posA); if (verbose) show_float_array(stderr, "A", &A, 1);
//     
//     free(A.e);
//   }
// 
// void test_float_array_sort_entries(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_sort_entries} ...\n");
//     
//     float_array_t A = float_array_new(11,17,0); 
//     float_array_throw(&A, 0.25);
//     float_array_scramble_entries(&A);
//     if (verbose) show_float_array(stderr, "A", &A, 10);
//     
//     auto void do_test(int orow, int ocol);
//     
//     void do_test(int orow, int ocol)
//       { if (verbose) fprintf(stderr, "ordering: orow = %+2d ocol = %+2d\n", orow, ocol);
//         float_array_scramble_entries(&A); 
//         float_array_sort_entries(&A, orow, ocol);
//         if (verbose) show_float_array(stderr, "A", &A, 10); 
//         check_order_of_entries(&A, orow, ocol);
//       }
//       
//     do_test(+2, +1);
//     do_test(+1, -2);
//     do_test(-1, 00);
//     free(A.e);
//   }
//   
// void test_float_array_copy(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_copy} ...\n");
//     
//     float_array_t A = float_array_new(11,17,0); 
//     float_array_throw(&A, 0.25); float_array_scramble_entries(&A);
//     
//     float_array_t R = float_array_new(25,200,0);
//     
//     float_array_copy(&A, &R);
//     check_float_array(&A, "A", &R, "R", verbose);
//     
//     free(A.e);
//     free(R.e);
//   }
//   
// void test_float_array_sort_entries_ins(int it, bool_t verbose)
//   {
//     float_array_t A = float_array_new(11,17,0); 
//     float_array_throw(&A, 0.25); 
//     float_array_scramble_entries(&A);
//     
//     
//     if (verbose) fprintf(stderr, "Checking {float_array_sort_entries_ins} ...\n");
//     float_array_t R = float_array_new(17,312,11);     
//     float_array_copy(&A, &R);
//     
//     float_array_sort_entries_ins(&A, +1, +2, 0,A.ents); if (verbose) show_float_array(stderr, "A", &A, 10);
//     float_array_sort_entries(&R, +1, +2); if (verbose) show_float_array(stderr, "R", &R, 10);
//     check_float_array(&A, "A", &R, "R", verbose);
//     
//     float_array_sort_entries_ins(&A, +2, +1, 0,A.ents); if (verbose) show_float_array(stderr, "A", &A, 10);
//     float_array_sort_entries(&R, +2, +1); if (verbose) show_float_array(stderr, "R", &R, 10);
//     check_float_array(&A, "A", &R, "R", verbose);
//   }
//   
// void test_float_array_condense(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_condense} ...\n");
//     
//     float_array_t A = float_array_new(0,0,0);
//     float_array_t R = float_array_new(0,0,0);
//     
//     int k, r;
//     int rMax = 5; /* Max repetitions of each index pair. */
//     float_array_pos_t posA = 0, posR = 0;
//     for (r = 0; r < rMax; r++)
//       { /* Generate some set of distinct pairs {i,j}: */
//         for (k = 0; k < 300; k++)
//           { int i = (k + (k*k/3)) / 19;
//             int j = (k + (k*k/3)) % 19;
//             int rCount = ((i + 3*j) % rMax) + 1; /* Number of repetitions for this {i,j}. */
//             if (r < rCount)
//               { double Aij = r;
//                 posA = float_array_add_elem(&A, posA, i, j, Aij);
//               }
//             if (r == 0)
//               { double Rij = rCount*(rCount-1)/2;
//                 posR = float_array_add_elem(&R, posR, i, j, Rij);
//               }
//           }
//       }
//     float_array_trim(&A, posA); float_array_sort_entries(&A, +2, +1); if (verbose) show_float_array(stderr, "A", &A, 10);
//     float_array_trim(&R, posR); float_array_sort_entries(&R, +2, +1); if (verbose) show_float_array(stderr, "R", &R, 10);
// 
//     auto double dbl_add(float_array_index_t row, float_array_index_t col, double v0, double v1);
//     
//     double dbl_add(float_array_index_t row, float_array_index_t col, double v0, double v1)
//       { return v0+v1; }
//     
//     float_array_condense(&A, dbl_add);
//     check_float_array(&A, "A", &R, "R", verbose);
//     
//     free(A.e);
//     free(R.e);
//   }
// 
// void test_float_array_mix(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_mix} ...\n");
//     
//     float_array_t A = float_array_new(0,0,0);
//     float_array_t B = float_array_new(0,0,0);
//     float_array_t R = float_array_new(0,0,0);
//     
//     int k;
//     float_array_pos_t posA = 0, posB = 0, posR = 0;
//     for (k = 0; k < 100; k++)
//       { int i = (k + (k*k/3)) / 19;
//         int j = (k + (k*k/3)) % 19;
//         double Aij = 1 + sin(i)*sin(3*j);
//         double Bij = 1 + cos(i)*cos(5*j);
//         double Rij = 0.5 * Aij + 2.0 * Bij;
//         posA = float_array_add_elem(&A, posA, i, j, Aij);
//         posB = float_array_add_elem(&B, posB, i, j, Bij);
//         posR = float_array_add_elem(&R, posR, i, j, Rij);
//       }
//     float_array_trim(&A, posA); float_array_sort_entries(&A, +2, +1); if (verbose) show_float_array(stderr, "A", &A, 10);
//     float_array_trim(&B, posB); float_array_sort_entries(&B, +2, +1); if (verbose) show_float_array(stderr, "B", &B, 10);
//     float_array_trim(&R, posR); float_array_sort_entries(&R, +2, +1); if (verbose) show_float_array(stderr, "R", &R, 10);
//     
//     float_array_t C = float_array_new(0,0,0);
//     float_array_mix(0.5, &A, 2.0, &B, &C);
//     check_float_array(&C, "C", &R, "R", verbose);
//     
//     free(A.e);
//     free(B.e);
//     free(C.e);
//     free(R.e);
//   }
//   
// void test_float_array_transpose(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_transpose} ...\n");
//     
//     float_array_t A = float_array_new(5,7,0); 
//     float_array_throw(&A, 0.25); 
//     float_array_scramble_entries(&A);
//     
//     float_array_transpose(&A);  if (verbose) show_float_array(stderr, "A", &A, 10);
//     free(A.e);
//   }
//   
// void test_float_array_extract_row_float_array_extract_col_float_array_mul(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array,float_array_extract_col,float_array_mul} ...\n");
//     
//     float_array_t A = float_array_new(7,11,0); 
//     float_array_throw(&A, 0.50); 
//     if (verbose) show_float_array(stderr, "A", &A, 10);
//     
//     float_array_t B = float_array_new(11,5,0); 
//     float_array_throw(&B, 0.50); 
//     if (verbose) show_float_array(stderr, "B", &B, 10);
//     
//     float_array_t C = float_array_new(0,0,0);
//     /* Multiply by the library, result is {C}: */
//     float_array_mul(&A, &B, &C); if (verbose) show_float_array(stderr, "C", &C, 10);
//     
//     /* Multiply by extracting rows and cols, result is {R}: */
//     double Arow[A.cols];
//     double Bcol[B.rows];
//     float_array_sort_entries(&A, +2, +1) /* (by row,col) */; if (verbose) show_float_array(stderr, "A", &A, 10);
//     float_array_sort_entries(&B, +1, +2) /* (by col,row) */; if (verbose) show_float_array(stderr, "B", &B, 10);
//     
//     float_array_t R = float_array_new(0,0,0);
//     float_array_pos_t posR = 0;
//     int i;
//     float_array_pos_t posA = 0;
//     for (i = 0; i < A.rows; i++)
//       { posA = float_array_extract_row(&A, posA, i, Arow, A.cols);
//         int j;
//         float_array_pos_t posB = 0; 
//         for (j = 0; j < B.cols; j++)
//           { posB = float_array_extract_col(&B, posB, j, Bcol, B.rows);
//             double sum = 0;
//             int k;
//             for (k = 0; k < A.cols; k++) { sum += Arow[k]*Bcol[k]; }
//             posR = float_array_add_elem(&R, posR, i,j, sum);
//           }
//       }
//     float_array_trim(&R, posR);
//     
//     check_float_array(&C, "C", &R, "R", verbose);
//     free(A.e);
//     free(B.e);
//     free(C.e);
//     free(R.e);
//   }
//   
// void test_float_array_map_row(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_map_row} ...\n");
//     
//     float_array_t A = float_array_new(7,11,0);
//     float_array_throw(&A, 0.50);
//     if (verbose) show_float_array(stderr, "A", &A, 10);
//     
//     float_array_dim_t nu = A.rows;
//     double u[nu];
//     float_array_dim_t nv = A.cols;
//     double v[nv];
//     int i,j;
//     for (i = 0; i < nu; i++) { u[i] = drandom(); }
//     /* By {float_array_map_row}, result in {v}: */
//     float_array_map_row(u, nu, &A, v, nv);
//     /* By the book, reult in {r}: */
//     float_array_dim_t nr = A.cols;
//     double r[nr];
//     float_array_sort_entries_ins(&A, +1, +2, 0, A.ents) /* (by col,row) */; if (verbose) show_float_array(stderr, "A", &A, 10);
//     double Acol[A.rows];
//     float_array_pos_t posA = 0;
//     for (j = 0; j < A.cols; j++) 
//       { double sum = 0; 
//         posA = float_array_extract_col(&A, posA, j, Acol, A.rows);
//         assert(nu == A.rows);
//         for (i = 0; i < nu; i++) { sum += u[i]*Acol[i]; }
//         r[j]= sum;
//       }
//     assert(posA == A.ents);
//     compare_vectors(v, "v", r, "r", nv, 1.0e-10);
//     free(A.e);
//   }
//   
// void test_float_array_map_col(int it, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_map_col} ...\n");
//     
//     float_array_t A = float_array_new(7,11,0);
//     float_array_throw(&A, 0.50);
//     if (verbose) show_float_array(stderr, "A", &A, 10);
//     
//     float_array_dim_t nu = A.cols;
//     double u[nu];
//     float_array_dim_t nv = A.rows;
//     double v[nv];
//     int i,j;
//     for (i = 0; i < nu; i++) { u[i] = drandom(); }
//     /* By {float_array_map_col}, result in {v}: */
//     float_array_map_col(&A, u, nu, v, nv);
//     /* By the book, reult in {r}: */
//     float_array_dim_t nr = A.rows;
//     double r[nr];
//     float_array_sort_entries(&A, +2, +1) /* (by row, col) */; if (verbose) show_float_array(stderr, "A", &A, 10);
//     double Arow[A.cols];
//     float_array_pos_t posA = 0;
//     for (i = 0; i < A.rows; i++) 
//       { double sum = 0; 
//         posA = float_array_extract_row(&A, posA, i, Arow, A.cols);
//         assert(A.cols == nu);
//         for (j = 0; j < nu; j++) { sum += Arow[j]*u[j]; }
//         r[i]= sum;
//       }
//     compare_vectors(v, "v", r, "r", nv, 1.0e-10);
//     free(A.e);    
//   }
//   
// void test_float_array_add_diagonal_float_array_fill_diagonal(int it, bool_t fill, bool_t verbose)
//   {
//     if (verbose) fprintf(stderr, "Checking {float_array_add_diagonal,float_array_fill_diagonal} ...\n");
//     
//     float_array_t C = float_array_new(17,31,0);
//     float_array_count_t nd = imin(lcm(C.cols, C.rows), 300);
//     double d[nd];
//     int k;
//     double vfill = drandom();
//     for (k = 0; k < nd; k++) { d[k] = (fill ? vfill : drandom()); }
//     float_array_pos_t posC = 0;
//     float_array_index_t row1 = 5, col1 = 3;
//     if (fill)
//       { posC = float_array_fill_diagonal(&C, posC, row1, col1, vfill, nd); }
//     else
//       { posC = float_array_add_diagonal(&C, posC, row1, col1, d, nd); }
//     float_array_trim(&C, posC);
//     
//     /* Check elements: */
//     float_array_t R = float_array_new(C.rows,C.cols,0);
//     float_array_pos_t posR = 0;
//     for (k = 0; k < nd; k++) 
//       { float_array_index_t row = (row1 + k) % R.rows;
//         float_array_index_t col = (col1 + k) % R.cols;
//         posR = float_array_add_elem(&R, posR, row, col, d[k]);
//       }
//     float_array_trim(&R, posR);
//          
//     check_float_array(&C, "C", &R, "R", verbose);
//     free(C.e);   
//   }
// 
// void float_array_scramble_entries(float_array_t *A)
//   {
//     float_array_pos_t p;
//     for (p = 0; p < A->ents; p++)
//       { float_array_pos_t q = p + (int)floor(drandom()*(A->ents - p));
//         if (q >= A->ents) { q = A->ents - 1; }
//         float_array_entry_t e = A->e[p];
//         A->e[p] = A->e[q];
//         A->e[q] = e;
//       }
//   }
  
// void check_order_of_entries(float_array_t *A, int orow, int ocol)
//   {
//     float_array_pos_t p1;
//     for (p1 = 1; p1 < A->ents; p1++)
//       { float_array_pos_t p0 = p1 - 1;
//         float_array_entry_t *e0 = &(A->e[p0]);
//         float_array_entry_t *e1 = &(A->e[p1]);
//         assert((e0->col < A->cols) && (e0->row < A->rows));
//         assert((e1->col < A->cols) && (e1->row < A->rows));
//         int cmp = array_compare_indices(e0->row, e0->col, e1->row, e1->col, orow, ocol);
//         if (cmp > 0) 
//           { fprintf(stderr, "** entries out of order (args = %+d %+d)\n", orow, ocol);
//             fprintf(stderr, "  A->e[%d] = [%d][%d] (%24.16e)\n", p0, e0->row, e0->col, e0->val);
//             fprintf(stderr, "  A->e[%d] = [%d][%d] (%24.16e)\n", p1, e1->row, e1->col, e1->val);
//             assert(FALSE);
//           }
//       }
//   } 



