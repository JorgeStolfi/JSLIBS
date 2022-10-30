#ifndef wt_table_pair_H
#define wt_table_pair_H

/* Weight tables for filtering digital signals */
/* Last edited on 2022-10-30 19:30:13 by stolfi */

#define wt_table_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <vec.h>
#include <bool.h>

/* MATCHED TABLES FOR MULTISCALE FILTERING

  The procedures in this section create two matched filter weight
  tables {*wtb0,*wtb1} for multiscale filtering. The first stage
  should use {wtb0}, and subsequent stages should use {wtb1}. They
  also put in {*wname0,*wname1} strings that describe those two
  tables.

  The tables will have odd length and will be normalized to unit sum.
  They are packaged as newly allocated {double_vec_t}s.
  
  !!! These procedures should be generalized to even lengths. !!! */

void wt_table_pair_make_gaussian
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  );
  /* Creates matched filter tables from Gaussian distributions. */
  
void wt_table_pair_make_binomial
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  );
  /* Creates matched filter tables from binomial distributions. */

void wt_table_pair_make_triangular
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  );
  /* Creates matched filter tables from triangular distributions. */

#endif
