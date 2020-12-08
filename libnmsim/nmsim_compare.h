#ifndef nmsim_compare_H
#define nmsim_compare_H
 
/* Value comparison for write/read tests. */
/* Last edited on 2019-07-04 15:58:32 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

void nmsim_compare_int64_param(char *name, int64_t vr, int64_t ve);
  /* Compares a parameter value {vr} read from a file
    with the expected value {ve}. Aborts if different. */
    
void nmsim_compare_double_param
  ( char *name, 
    double vr, 
    double ve, 
    double prec,
    bool_t special_0, 
    bool_t special_1
  );
  /* Compares a parameter value {vr} read from a file with the
    expected value {ve}. Requires identity if either is {NAN} or
    infinite. 
    
    If {special_0} is true, requires identity if either of the two
    values is 0.0. If neither is zero, but {abs(ve)} is less than
    {2*prec}, requires {abs(vr)} less than {2*prec} and same sign.
    
    If {special_1} is true, requires identity if either of them is
    {+1.0} or {-1.0}.  If neither is {+1.0} or {-1.0},
    and {abs(1-abs(ve))} is less than {2*prec},
    requires {abs(abs(vr)-abs(ve))} less than {2*prec}
    
    Otherwise, aborts if abs difference is more than
    {prec}. */
    
#endif
