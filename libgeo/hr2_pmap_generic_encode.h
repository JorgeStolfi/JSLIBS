#ifndef hr2_pmap_generic_encode_H
#define hr2_pmap_generic_encode_H

/* Last edited on 2024-09-17 16:29:39 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>
#include <hr2_pmap.h>
#include <sign.h>

void hr2_pmap_generic_encode(hr2_pmap_t *M, double y[]);
  /* Converts a generic projective map {M} to nine parameters {y[0..8]},
    which are the elements of {M.dir} linearized by rows,
    except that {y[3..5]} and {y[6..8]} are swapped if the determinant of {M.dir}
    is negative.  */

void hr2_pmap_generic_decode(double y[], hr2_pmap_t *M);
  /* Sets the elements of {M.dir}, in row-by-row order,
    to {y[0..8]}.  Then swaps the last two rows, if necessary
    to make the handedness of the map be positive.
    Then computes {M.inv}.  */

#endif
