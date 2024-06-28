#define PROG_NAME "test_jsdebug"
#define PROG_DESC "test of {jsmath.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-06-28 02:17:47 by stolfi */ 
/* Created on 2011-09-20 by J. Stolfi, UNICAMP */

#define test_jsdebug_COPYRIGHT \
  "Copyright © 2011  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsmath.h>
#include <jsfile.h>
#include <affirm.h>
#include <jsstring.h>
#include <bool.h>

#include <jsdebug.h>

int32_t main (int32_t argn, char **argv)
  { 
    fprintf(stderr, "  sizeof(string_t) = %lu\n", sizeof(string_t));
    
    int32_t nm = 50;
    string_t mat[nm];
    
    jsdebug_addr_span("mat", &(mat[0]), &(mat[nm]), nm);

    return 0;
  }
