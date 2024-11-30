/* test_gausol_print.c --- test program for gausol_printularize.h  */
/* Last edited on 2024-11-29 07:18:36 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <jsmath.h>
#include <jsrandom.h>
#include <bool.h>
#include <affirm.h>

#include <gausol_print.h>
#include <gausol_test_tools.h>

int32_t main(int32_t argc, char **argv);

void test_gausol_print_row_has_name_eq(bool_t verbose);
void test_gausol_print(uint32_t trial, bool_t verbose);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {

    test_gausol_print_row_has_name_eq(FALSE);
    uint32_t ntrials = 20;

    for (uint32_t trial = 0;  trial < ntrials; trial++) 
      { test_gausol_print(trial, trial < 10); }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }
  
void test_gausol_print_row_has_name_eq(bool_t verbose)
  {
    for (uint32_t m = 5; m <= 6; m++)
      { for (uint32_t rank = 0; rank <= m; rank++)
          { fprintf(stderr, "\n");
            fprintf(stderr, "m = %d rank= %d\n\n", m, rank);
            for (uint32_t i = 0; i < m; i++)
             { bool_t at_dash, at_data;
                gausol_print_row_has_name_eq(i, m, rank, &at_dash, &at_data);
                if ((rank >= 1) && (rank <= m-1) && (i == rank))
                  { fprintf (stderr, " %s ---\n", (at_dash ? "NAME =" : "      ")); }
                fprintf (stderr, " %s %03d\n", (at_data ? "NAME =" : "      "), i);
              }
            fprintf(stderr, "\n");
          }
      }
    fprintf(stderr, "done.\n");
  }
       

void test_gausol_print (uint32_t trial, bool_t verbose)
  {
    srand(20230225 + 2*trial);
    srandom(20230225 + 2*trial);

    uint32_t m = (trial+3/4);
    uint32_t n = uint32_abrandom(3, 5);
    uint32_t p = uint32_abrandom(0, 3);
    uint32_t q = uint32_abrandom(0, 3);

    uint32_t rank = uint32_abrandom(0, (m < n ? m : n));
    
    uint32_t *prow = (drandom() < 0.25 ? NULL : random_perm(m, NULL));
    uint32_t *pcol = (drandom() < 0.25 ? NULL : random_perm(n, NULL));

    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "%s  trial = %d\n", __FUNCTION__, trial);
    fprintf(stderr, "testing with m = %d  n = %d  p = %d  q = %d", m, n, p, q);
    fprintf(stderr, " rank = %d ...\n", rank);
    
    double tiny = 1.0e-2;
    double pzero = 0.25;
    double mag = 99.999;
    double *A = gausol_test_tools_throw_matrix(m, n, pzero, mag, mag, tiny, "", "A", verbose);
    double *B = gausol_test_tools_throw_matrix(m, p, pzero, mag, mag, tiny, "", "B", verbose);
    double *C = gausol_test_tools_throw_matrix(m, q, pzero, mag, mag, tiny, "", "C", verbose);
    
    fprintf(stderr, "--- tests without permutations ---\n");
    gausol_print_system(stderr, 6, "%8.3f", "matrices A,B,C",    m, NULL,rank,  n, NULL,rank,  "A",   A,  p,  "B",   B,  q, "C",C, "footer");
    gausol_print_system(stderr, 6, "%8.3f", "matrices A,C only", m, NULL,rank,  n, NULL,rank,  "A",   A,  0, NULL,NULL,  q, "C",C, "footer");
    gausol_print_system(stderr, 6, "%8.3f", "matrices B,C only", m, NULL,rank,  0, NULL,0, NULL,NULL,  p,  "B",   B,  q, "C",C, "footer");
    
    fprintf(stderr, "--- tests with row perms only ---\n");
    gausol_print_system(stderr, 6, "%8.3f", "matrices A,B,C",    m, prow,rank,  n, NULL,rank,  "A",   A,  p,  "B",   B,  q, "C",C, "footer");
    
    fprintf(stderr, "--- tests with column perms only ---\n");
    gausol_print_system(stderr, 6, "%8.3f", "matrices A,B,C",    m, NULL,rank,  n, pcol,rank,  "A",   A,  p,  "B",   B,  q, "C",C, "footer");
    
    fprintf(stderr, "--- tests with row and column perms  ---\n");
    gausol_print_system(stderr, 6, "%8.3f", "matrices A,B,C",    m, prow,rank,  n, pcol,rank,  "A",   A,  p,  "B",   B,  q, "C",C, "footer");
    
    fprintf(stderr, "done.\n");
    fprintf(stderr, "======================================================================\n");
  }
