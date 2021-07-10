#ifndef wt_table_H
#define wt_table_H

/* Weight tables for filtering digital signals */
/* Last edited on 2021-07-04 09:36:03 by jstolfi */

#define wt_table_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

/* WEIGHT TABLES WITHOUT ALLOCATION

  The procedures in this section create filter weight tables of 
  arbitrary length {n>0}, even or odd, stored into client-given arrays.  */

void wt_table_fill_gaussian(double sigma, int n, double wt[]);
  /* Stores into {wt[0..n-1]} a table derived from a Gaussian
    distribution {G(z)} with mean {n/2} and standard deviation {sigma}.
    
    More precisely, {wt[i]} will be the integral of {G(z)} over 
    the interval {[i _ i+1]}. The table is normalized so that
    the sum of all entries is 1. */

void wt_table_fill_binomial(int n, double wt[]);
  /* Sets {wt[i] = choose(i,n-1)/2^(n-1)} for {i} in {0..n-1}. 
    Note that the sum of all entries is 1.
    The variance will be {(n-1)/4}, so the
    standard deviation will be {sqrt(n-1)/2)}.
    
    The table will be a partition of unity if used with stride
    1 or (if n >=2) with stride 2. */

void wt_table_fill_triangular(int n, double wt[]);
  /* First, builds a triangular weight table, {wt[i] = 1 - |i-c|/h}
    where {c = (n-1)/2} and {h=(n+1)/2}for {i} in {0..n-1}. 
    Requires {n} odd. 
    
    Thus, if {n} is odd, the central element will be set to 1, and the
    table is the convolution of two uniform distributions on {r+1}
    elements. In any case, the values will be symmetric; the elements
    {wt[-1]} and {wt[n]}, if they existed, would be zero, and {wt[0] =
    wt[n-1] = 2/(n+1)}.
    
    The table is then normalized so that the sum of all entries is 1. If
    {n} is odd, {2*r+1}, the table will be a partition of unity when
    used with stride {r+1}. */

void wt_table_fill_hann(int n, double wt[]);
  /* Stores into {wt[0..n-1]} a Hann-like weight table.
  
    First, sets each {wt[i]} to {0.5*(1 + cos(pi*(i-c)/h))} where
    {c=(n-1)/2} and {h=(n+1)/2}. Thus, if {n} is odd, the central element
    will be set to 1. In any case, the values will be symmetric; the
    elements {wt[-1]} and {wt[n]}, if they existed, would be zero, and
    {wt[0] = wt[n-1] = 0.5*(1 + cos(pi*(n-1)/(n+1))}.
    
    The table is then normalized so that the sum of all entries is 1.
    If {n} is odd, {2*r+1}, the table will be a partition of unity 
    when used with stride {r+1}. */

/* WEIGHT TABLES WITH ALLOCATION

  The procedures in this section create weight tables of 
  odd length {n = 2*r+1}, packaged as newly allocated {double_vec_t}s. 
  
  !!! These procedure should be generalized to even lengths. !!! */

double_vec_t wt_table_make_gaussian(double sigma, double maxLoss);
  /* Builds a weight table derived from a Gaussian distribution {G(z)}
    with standard deviation {sigma}.
    
    The procedure chooses {r} so that the total omitted (out-of-window)
    probability (the integral of {G(z)} outside the interval
    {[-(r+1/2) _ +(r+1/2)]}) does not exceed {maxLoss}.
    The table will have length {n=2*r+1}, and will be filled
    as in {wt_table_fill_gaussian}. */

double_vec_t wt_table_make_binomial(int r);
  /* Builds a weight table with length {n=2*r+1}, derived from a
    binomial distribution, as in {wt_table_fill_binomial}. The
    variance will be {r/2}, so the standard deviation will be
    {sqrt(r/2)}. */

double_vec_t wt_table_make_triangular(int r);
  /* Builds a weight table with length {n=2*r+1} and a
    triangular value profile, as in {wt_table_fill_triangular}.
    !!! What is the variance? !!! */

double_vec_t wt_table_make_hann(int r);
  /* Builds a weight table with size {n=2*r+1} and a
    Hann value profile, as in {wt_table_fill_hann}.
    !!! What is the variance? !!!. */

/* WEIGHT TABLE ATTRIBUTES */

double wt_table_avg(int n, double wt[]);
  /* Computes the average index of{wt[0..n-1]}, that is,
    {SUM{i*wt[i]}/SUM{wt[i]}}.  Assumes {wt[i]} is non-negative;
    Returns {NaN} if the sum of {wt} is zero. */

double wt_table_var(int n, double wt[], double avg);
  /* Computes the variance of the indices of{wt[0..n-1]}
    with respect to the assumed average index {avg},
    that is, {SUM{(i-avg)^2*wt[i]}/SUM{wt[i]}}.  
    Assumes {wt[i]} is non-negative;  Returns {NaN}
    if the sum of {wt} is zero. */

/* PRINTOUT */

void wt_table_print(FILE *wr, char *wtname, int n, double wt[], int stride);
  /* prints to {wr} the name {wtname} (if not NULL), the weight table
    {wt[0..n-1]}, and its main statistical properties. 
    
    If {stride} is positive, also prints the effect of overlapping 
    copies of the table with displacement {stride}. */

char *wt_table_make_descr(int n, double wt[], char *fmt);
  /* Returns the weight table {wt[0..n-1]} formatted as a character string,
    with each element printed with the given {fmt}. */

/* PARSING COMMAND LINE ARGUMENTS */

double_vec_t wt_table_args_parse(argparser_t *pp, bool_t unitNorm);
   /* Parses a sequence of numbers from the command line, in the
     format described by {wt_table_args_HELP} below. \
     
     If {unitNorm} is TRUE, scales all weights so that their sum is 1;
     the optional "/ {DENOM}" is parsed but ignored. 
     
     If {unitNorm} is FALSE, does no normalization; but, if the
     optional "/ {DENOM}" is present, divides all weights by
     {DENOM}. */
     
#define wt_table_args_HELP \
  "{WEIGHT}.. [ \"/\" {DENOM} ]"

#define wt_table_args_INFO \
  "The weights may be integer or fractional," \
  " and may include \"e\" exponent parts. If the" \
  " optional slash and denominator are present," \
  " all weights are divided by {DENOM}."

#define wt_table_args_norm_sum_INFO \
  "The weights may be integer or fractional," \
  " and may include \"e\" exponent parts. They are" \
  " implicitly normalized to unit sum. (The slash" \
  " and {DENOM}, if present, are parsed but ignored.)"

/* INTERNAL PROTOTYPES */

bool_t wt_table_check_normalization(int n, double wt[], double tol, bool_t die);
  /* Checks whether {wtb.e[0..wtb.ne-1]} add to 1. If not, either
    dies with error (if {die} is TRUE) or returns FALSE silently 
    (if {die} is FALSE). */

double wt_table_gaussian_loss(int n, double sigma);
  /* Returns the weight that is lost if a table with {n} entries is
    used to represent a Gaussian weight distribution with standard
    deviation {sigma}. */

double wt_table_gaussian_entry(int n, int k, double sigma);
  /* Computes the entry with position {k} in a table with {n}
    elements, where each element is the integral of a truncated
    Gaussian weight distribution with standard deviation {sigma}
    in a unit sub-interval of the interval {[-n/2 _ +n/2]}. The
    position {k} ranges in {0..n-1}. The table is normalized
    to unit sum. 
    
    !!! Inefficient -- should reuse the shared {erf} calls !!! */

#endif
