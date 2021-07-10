#ifndef dnae_weight_H
#define dnae_weight_H

/* Weight tables for filtering digital signals */
/* Last edited on 2014-06-10 14:01:42 by stolfilocal */

/* !!! Superseded by {wt_table.h} in JSLIBS/libjs? !!! */

#define dnae_weight_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>

#include <vec.h>
#include <argparser.h>

void dnae_weight_table_print(FILE *wr, char *wtname, double_vec_t *wt);
  /* Prints out the weight table {wt} (ant its name {wtname}) 
    to file {wr}. */

char *dnae_weight_make_descr(double_vec_t *wt, char *fmt);
  /* Returns the weight table {wt} formatted as a character string,
    with each element printed with the given {fmt}. */

double_vec_t dnae_weight_table_make_gaussian(double sigma, double maxLoss);
  /* Builds a weight table derived from a Gaussian distribution {G(z)}
    of standard deviation {sigma}. The table length will be odd,
    {2*r+1} for some {r>=0}. The procedure chooses {r} so that the
    total omitted (out-of-window) probability does not exceed
    {maxLoss}.
    
    More precisely, {w[i]} will be the integral of {G(z)} within a
    unit interval centered on {i-r}, for {i = 0..2*r}. The omitted
    probability is the integral of {G(z)} outside the interval
    {[-(r+1/2) _ +(r+1/2)]}.
    
    The table is normalized so that the sum of all entries is 1. */

double_vec_t dnae_weight_table_make_binomial(int r);
  /* Builds a weight table derived from a binomial distribution {w[i]
    = choose(i,n)/2^n} for {i} in {0..n}, where {n = 2*r}. Note that
    this is the convolution of {2*r} binary distributions. The table
    length will be {2*r+1}. The variance will be {r/2}, so the
    standard deviation will be {sqrt(r/2)}. */

double_vec_t dnae_weight_table_make_triangular(int r);
  /* Builds a triangular weight table, {w[i] = (r+1 - |i-r|)/(r+1)^2}
    for {i} in {0..2*r}. Note that this is the convolution of two
    uniformdistributions on {r+1} elements. The table length will be
    {2*r+1}. */

/* MATCHED TABLES FOR MULTISCALE FILTEING

  The procedures in this section create two matched filter weight tables
  {*wtb0,*wtb1} for multiscale filtering. The first
  stage should use {wtb0}, and subsequent stages should use {wtb1}. 
  They also put in {*wname0,*wname1} strings that describe those
  two tables. */

void dnae_weight_table_pair_make_gaussian
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  );
  /* Creates matched filter tables from Gaussian distributions. */
  
void dnae_weight_table_pair_make_binomial
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  );
  /* Creates matched filter tables from binomial distributions. */

void dnae_weight_table_pair_make_triangular
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  );
  /* Creates matched filter tables from triangular distributions. */

/* PARSING COMMAND LINE ARGUMENTS */

double_vec_t dnae_weight_args_parse(argparser_t *pp, bool_t unitNorm);
   /* Parses a sequence of numbers from the command line, in the
     format described by {dnae_weight_args_HELP} below. \
     
     If {unitNorm} is TRUE, scales all weights so that their sum is 1;
     the optional "/ {DENOM}" is parsed but ignored. 
     
     If {unitNorm} is FALSE, does no normalization; but, if the
     optional "/ {DENOM}" is present, divides all weights by
     {DENOM}. */
     
#define dnae_weight_args_HELP \
  "{WEIGHT}.. [ \"/\" {DENOM} ]"

#define dnae_weight_args_INFO \
  "The weights may be integer or fractional," \
  " and may include \"e\" exponent parts. If the" \
  " optional slash and denominator are present," \
  " all weights are divided by {DENOM}."

#define dnae_weight_args_norm_sum_INFO \
  "The weights may be integer or fractional," \
  " and may include \"e\" exponent parts. They are" \
  " implicitly normalized to unit sum. (The slash" \
  " and {DENOM}, if present, are parsed but ignored.)"

#endif
