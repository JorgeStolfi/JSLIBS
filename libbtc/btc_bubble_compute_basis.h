#ifndef btc_bubble_compute_basis_H
#define btc_bubble_compute_basis_H

/* Computing the basis of a bubble model for the BTC price. */
/* Last edited on 2015-04-20 00:23:01 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_compute_basis(int nd, int nb, btc_bubble_t bp[], int hrad, double bval[]);
  /* Fills the array {bval} with the bubble function basis for 
    {nd} days and {nb} bubbles with parameters {bp[0..nb-1]}.

    Namely, sets {bval[nb*id + jb]} to the value of bubble function {jb}
    on day {id}, for {jb} in {0..nb-1} and {id} in {0..nb-1}. The bubble
    functions are smoothed with Hann window of radius {hrad}. Each
    bubble function has the maximum value 1.0 (before smoothing). */
  

#endif
