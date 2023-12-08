#ifndef wt_table_H
#define wt_table_H

/* Weight tables for filtering digital signals */
/* Last edited on 2023-11-25 18:48:44 by stolfi */

#define wt_table_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>

/* INDEX STATISTICS */

double wt_table_avg(int32_t n, double wt[]);
  /* Computes the weighted average INDEX of{wt[0..n-1]}, that is,
    {SUM{i*wt[i]}/SUM{wt[i]}}.  Assumes {wt[i]} is non-negative;
    Returns {NaN} if the sum of {wt} is zero. */

double wt_table_var(int32_t n, double wt[], double avg);
  /* Computes the variance of the indices of{wt[0..n-1]} with respect to
    the assumed average index {avg}, that is,
    {SUM{(i-avg)^2*wt[i]}/SUM{wt[i]}}. Assumes {wt[i]} is non-negative;
    Returns {NaN} if the sum of {wt} is zero. */

/* NORMALIZATION */

void wt_table_normalize_sum(int32_t n, double wt[]);
  /* Scales the elements of {wt[0..n-1]} so that their sum is 1. */

/* SHIFTED SUMMATION */

void wt_table_shifted_sum(int32_t n, double wt[], int32_t stride, double ws[]);
  /* Returns in {ws[0..n-1]} the result of adding all copies of {wt}
    shifted by integer multiples of {stride}. Namely, sets {ws[k]} to
    the sum of {wt[k+i*stride]} for all integer {i}, assuming that
    non-existent elements of {wf} are zero. Requires {stride >= 1}. */

double_vec_t wt_table_convolution(int32_t n1, double wt1[], int32_t n2, double wt2[], int32_t stride);
  /* Returns a weight table {ws} that is the result of adding copies of {wt1} 
    shifted by multiples of {stride} and scaled according to {wt2}.
    
    More precisely, the returned table will have {ns = n1 + n2*stride} elements,
    and, for each {i} in {0..ns-1, {ws[i]} will be
    {SUM{wt1[i - k*stride]*wt2[k] : k\in 0..n2-1 }}, assuming that
    non-existent elements of {wt1} are zero. Requires {stride >= 1}. */

/* PRINTOUT */

void wt_table_print(FILE *wr, char *wtname, int32_t n, double wt[], int32_t stride);
  /* prints to {wr} the name {wtname} (if not NULL), the weight table
    {wt[0..n-1]}, and its main statistical properties. 
    
    If {stride} is positive, also prints the effect of overlapping 
    copies of the table with displacement {stride}. */

char *wt_table_make_descr(int32_t n, double wt[], char *fmt);
  /* Returns the weight table {wt[0..n-1]} formatted as a character string,
    with each element printed with the given {fmt}. */
/* TESTING AND DEBUGGING */

bool_t wt_table_check_normalization(int32_t n, double wt[], double tol, bool_t die);
  /* Checks whether {wtb.e[0..wtb.ne-1]} add to 1, with tolerance {tol}.
    If not, either dies with error (if {die} is {TRUE}) or returns
    {FALSE} silently (if {die} is {FALSE}). */

bool_t wt_table_check_partition_of_constant
  ( int32_t n, 
    double wt[], 
    int32_t stride,
    double tol, 
    bool_t die
  );
  /* Checks whether the sum of copies of {wt} shifted by all integer
    multiples of {stride} (which must be positive) is a sequence with
    all elements equal, within the tolerance {tol}. If not,
    either dies with error (if {die} is {TRUE}) or returns {FALSE} silently (if
    {die} is {FALSE}). */

#endif
