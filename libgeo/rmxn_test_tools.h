/* rmxn_test_tools.h --- tools for testing the {libgeo} functions. */
/* Last edited on 2024-12-05 10:28:41 by stolfi */
/* Created 2005-07-20 by J. Stolfi. */

#ifndef rmxn_test_tools_H
#define rmxn_test_tools_H

#include <stdio.h>
#include <stdint.h>

#include <affirm.h>
#include <rn.h>

void rmxn_test_tools_check_all_different(uint32_t m, uint32_t n, double *A, char *msg);
  /* Compares all pairs of elements of {A[0..m*n-1]} and bombs out if
    two are equal.  Used to check arrays randomly sampled
    from a distribution where any such identity is extremely unlikely,
    and most probably indicates an error in the generator. */

#endif
