/* See {multifok_score.h}. */
/* Last edited on 2023-02-12 06:52:35 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <wt_table.h>
#include <affirm.h>
#include <fget.h>
#include <bool.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>

#include <multifok_score.h>

double multifok_score_from_terms(int32_t NT, double wt[], double term[])
  { double score = 0.0;
    for (int32_t kt = 0; kt < NT; kt++)
      { /* Compute and combine term values: */
        score += wt[kt]*term[kt];
      }
    return score;
  }

void multifok_score_read_term_names_and_weights
  ( FILE *rd, 
    int32_t NB, 
    char *belName[],
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P,
    int32_t *NT_P, 
    double **wt_P, 
    char ***termName_P,
    bool_t verbose
  )
  {
    /* Read term names and weights from file: */
    int32_t NT = 0;
    string_vec_t termName = string_vec_new(50);
    double_vec_t wt = double_vec_new(50);
    while (TRUE)
      { bool_t ok = fget_test_comment_or_eol(rd, '#');
        if (ok) { continue; }
        if (fget_test_eof(rd)) { break; }
        /* There is something there: */
        char *tnk = fget_string(rd);
        double wtk = fget_double(rd);
        fget_comment_or_eol(rd, '#');

        /* Save in tables: */
        int32_t kt = NT;
        if (verbose) { fprintf(stderr, "%3d  %+16.12f %s", kt, wtk, tnk); }
        string_vec_expand(&termName,NT);
        double_vec_expand(&wt,NT);
        termName.e[kt] = tnk;
        wt.e[kt] = wtk;
        NT++;
      }
    string_vec_trim(&termName,NT);
    double_vec_trim(&wt,NT);
    
    /* Convert term names to indices: */
    int32_t NP;
    multifok_term_prod_t *prix;
    multifok_term_indices_from_names(NB, belName, NT, termName.e, &NP, &prix, verbose);
    
    if (verbose)
      { fprintf(stderr, "product index table:\n");
        for (int32_t ip = 0; ip < NP; ip++)
          { multifok_term_prod_t *pri = &(prix[ip]);
            fprintf(stderr, "%3d  %3d %3d  %3d %s\n", ip, pri->jb1, pri->jb2, pri->kt, pri->name);
          }
      }
    
    (*NT_P) = NT;
    (*termName_P) = termName.e;
    (*prix_P) = prix;
    (*wt_P) = wt.e;
  }

void multifok_score_write_term_names_and_weights(FILE *wr, int32_t NT, char *termName[], double wt[])  
  { for (int32_t kt = 0; kt < NT; kt++)
      { fprintf(wr, "%3d %+12.4f %s\n", kt, wt[kt], termName[kt]); }
    fflush(wr);
  }

#define multifok_score_C_COPYRIGHT \
    "© 2022 by the State University of Campinas (UNICAMP)"

