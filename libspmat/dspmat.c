/* See {dspmat.h}. */

#define dspmat_C_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2009-01-17 18:59:41 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <bool.h>

#include <spmat.h>
#include <spmat_io.h>
#include <spmat_linalg.h>

#include <dspmat.h>

#define dspmat_trivial_elem (0.0)
#define dspmat_elem_is_trivial(X) ((X)==0.0)

spmat_impl(dspmat_t, dspmat, double);

void dspmat_elem_write(FILE *wr, double *valP) { fprintf(wr, "%24.16e", *valP); }
void dspmat_elem_read(FILE *rd, double *valP) { (*valP) = fget_double(rd); }

spmat_io_impl(dspmat_t, dspmat, double);

#define dspmat_elem_add(X,Y) ((X)+(Y))
#define dspmat_elem_mul(X,Y) ((X)*(Y))
#define dspmat_elem_one (1.0)

spmat_linalg_impl(dspmat_t, dspmat, double);
