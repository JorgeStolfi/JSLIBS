#ifndef wt_table_pair_H
#define wt_table_pair_H

/* Weight tables for filtering digital signals */
/* Last edited on 2022-10-31 14:11:39 by stolfi */

#define wt_table_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <vec.h>
#include <bool.h>

/* MATCHED TABLES FOR MULTISCALE FILTERING

  The procedures in this section create two matched filter weight tables
  {wt0,wt1} for multiscale filtering. The first stage should use
  {wt0}, and subsequent stages should use {wt1}. These tables are
  returned in {*wt0_P,*wt1_P}. The procedures also put in
  {*wname0_P,*wname1_P} strings that describe those two tables.

  The tables will have odd length and will be normalized to unit sum.
  Their variances will be at least {var0} and {var1}, respectively.
  They are packaged as newly allocated {double_vec_t}s.
  
  If {verbose} is true, each procedure prints the two tables and 
  their startistical properties.
  
  !!! These procedures should be generalized to even lengths. !!! */

void wt_table_pair_make_gaussian
  ( double var0,
    double_vec_t *wt0_P, 
    char **wname0_P, 
    double var1,
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  );
  /* Creates matched filter tables from Gaussian distributions. */
  
void wt_table_pair_make_binomial
  ( double var0,
    double_vec_t *wt0_P, 
    char **wname0_P, 
    double var1,
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  );
  /* Creates matched filter tables from binomial distributions. */

void wt_table_pair_make_triangular
  ( double var0,
    double_vec_t *wt0_P, 
    char **wname0_P, 
    double var1,
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  );
  /* Creates matched filter tables from triangular distributions. */

#endif
