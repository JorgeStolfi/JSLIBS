/* Last edited on 2025-04-27 02:57:30 by stolfi */
/* Computing deltas and weights for the star-cycle transformation. */ 

#ifndef psgr_shrink_stacy_H
#define psgr_shrink_stacy_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <psgr_types.h>
#include <psgr.h>

void psgr_shrink_stacy_get_star_data
  ( psgr_t *gr,
    psgr_arc_t a,
    uint32_t n,
    double d_star[],
    double w_star[],
    double *wtot_P
  );
  /* Enumerates the arcs in the star (onext ring) of the origin vertex {org} of {a},
    which should have outdegree {n}.  Saves their deltas and weights in {d_star[0..n-1]}
    and {w_star[0..n-1]}, respectively.  Specifically, {d_star[k]} and {w_star[k]}
    will be the delta and weight of the arc {onext^k(a)}.
    
    Also computes the total weight {wtot} ad returns it in {*wtot_P}.
    
    The weights must be all strictly positive. */

void psgr_shrink_stacy_compute_cycle_data
  ( psgr_t *gr,
    uint32_t n,
    double d_star[],
    double w_star[],
    double wtot,
    double d_cycle[],
    double w_cycle[]
  );
  /* Given the weights and deltas {w_star[0..n-1],d_star[0..n-1]} of the
    arcs in the star of a vertex, which should have outdegree {n}, and
    the total weight {wtot}, computes the weights and deltas
    {w_cycle[0..n-1],d_cycle[0..n-1]} of the arcs in the
    quasi-equivalent cycle.
    
    Specifically, assuming {d_star[k]} and {w_star[k]} are the delta and
    weight of arc {an[k]=onext^k(a)} for some arc {a}, the procedure sets
    {d_cycle[k]} and {w_cycle[k]} to the delta and weight of the new
    edge to be added between {dst(an[k])} and {dst(an[k+1])} (where the arc indices
    are taken modulo {n}).
    
    The input star weights {w_star[0..n-1]} must be all positive, but some 
    resulting weights {w_cycle[0..n-1]} may be zero because of underflow. */
    
#endif
