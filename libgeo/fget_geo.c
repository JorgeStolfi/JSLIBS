/* See fget_geo.h. */
/* Last edited on 2021-06-09 19:48:33 by jstolfi */

/* Copyright © 2008 Jorge Stolfi, Unicamp. See note at end of file. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <r6.h>
#include <fget.h>
#include <affirm.h>

#include <fget_geo.h>

void fget_rn(FILE *rd, double p[], int32_t n)
  { int32_t i;
    for (i = 0; i < n; i++)
      { p[i] = fget_double(rd); }
  }

r2_t fget_r2(FILE *rd)
  { r2_t p;
    fget_rn(rd, p.c, 2);
    return p;
  } 

r3_t fget_r3(FILE *rd)
  { r3_t p;
    fget_rn(rd, p.c, 3);
    return p;
  } 

r4_t fget_r4(FILE *rd)
  { r4_t p;
    fget_rn(rd, p.c, 4);
    return p;
  } 
  
r6_t fget_r6(FILE *rd)
  { r6_t p;
    fget_rn(rd, p.c, 6);
    return p;
  } 
  
r3_t fget_r3_dir(FILE *rd)
  { r3_t d;
    fget_rn(rd, d.c, 6);
    r3_dir(&d, &d);
    return d;
  } 

/* Copyright © 2008 by Jorge Stolfi.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appears in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty of any kind.
*/
 
