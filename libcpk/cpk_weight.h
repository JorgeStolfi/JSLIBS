#ifndef cpk_weight_H
#define cpk_weight_H

/* Weight formula for evaluating proposed auction sets. */
/* Last edited on 2024-12-31 11:07:23 by stolfi */

#include <stdint.h>
#include <stdio.h>

/* 
  WEIGHT FORMULA
  
  Potential solutions are ranked according to their /total weight/,
  the sum { W(J) = W[J[i]] : i = 0..nJ-1 }. Here {W[v]} is the /weight/ or
  `desirablilty' of holding an auction centered at {v}, ignoring 
  the incompatibilty constraints.
  
  In this package, the weight {W[v]} is derived from four terms:
  
    {NUM(v)} - always 1. This term contributes {1*nJ} to the total
      weight {W(J)}, and thus favors solutions that have as many
      auction disks as possible.
      
    {URB(v)} - the expected audience of the winning station. More
      precisely, let {f(u)} be the fraction of the disk with center
      {u} and radius {rRad = dMin/2} that lies within the urbanized
      area; {URB(v)} is the expected value of {f(u)}, for any point
      within the auction disk centered at {v}. In particular, this
      term is 1 if the circle with radius {rAuc+rRad} centered at {v}
      is entirely within the urban area. This term favors solutions
      that maximize the fraction of the urban area that is serviced by
      CR stations.
      
    {DEM(v)} - the known demand for an auction centered at {v}. More
      precisely, the number of known DIs that would be eligible to
      enter an auction centered at {v}; that is, the number of DIs
      within distance {rAuc} from {v}. This term favors solutions that
      maximize the number of DIs that can compete for the auctions.
      
    {CLS(v)} - the nearness to known demand points.  More precisely,
      let {h} be the function {h(r) = max(0, 1-(r/rAuc)^2)};
      then {CLS(v)} is the sum of {h(dist(v,p))} over all 
      DIs {p}. This term favors solutions where the auction
      centers are as close to the DIs as possible.
      
  Note that each term ranges from 0 to the maximum number of DIs that
  can be covered by a disk or radius {rAuc}. The influence of each
  term on the total weight {W} is defined by a /priority/ {p} in the
  range {0..4} and a /coefficient/ {c} in the range {0 _ 1}. The
  coefficients {c}s of same priority are normalized so that their sum
  is 1.
  
  Terms with priority 4 or more are simply ignored. Terms with the
  same priority {p} in {0..3} are multiplied by their {c}s and added
  together to form a partial weight {Z[p]}. Then the sums {Z[0..3]}
  are combined into a single weight {W} in such a way that {Z[0]} will
  dominate the ranking of solutions; {Z[1]} will be used to break ties
  in {Z[0]}; and similarly for {Z[2]} and {Z[3]}.
  
  This last step is somewhat compromised by the limited precision of a
  double-precision number (roughly 15 decimal digits). Thus, small
  differences in {URB(v)} or {CLS(v)} may be lost to rounding.
  Moreover, when the number of auction points {nJ} and/or the number
  of DIs is very large (over 1000), there may be spill-over from one
  level to the next, so that a very large difference in {Z[p]} may
  overcome a small contrary difference in {Z[p-1]}; and {Z[3]} may be
  partly lost. */

typedef struct cpk_wcoeffs_t
  {
    uint32_t pNUM; double cNUM; /* For number of auctions. */
    uint32_t pURB; double cURB; /* For urbanization of service area. */
    uint32_t pDEM; double cDEM; /* For number of eligible DIs. */
    uint32_t pCLS; double cCLS; /* For closeness to DIs. */
  } cpk_wcoeffs_t;
  /* 
    Priorities and coefficients of terms in the weight formula. */

double cpk_compute_weight(cpk_wcoeffs_t *wc, double URB, double DEM, double CLS);
  /* Combines the terms {NUM} (assumed to be 1), {URB}, {DEM}, and
    {CLS} into a single wheight {W}, as specified by the priorities
    and coefficients in {wc}. */

void cpk_normalize_weight_coeffs(cpk_wcoeffs_t *wc);
  /* Scales the coefficients {wc->cNUM}, {wc->cURB}, etc. so that the
    coeffs with same priority {p} in {0..3} add to 1 (maintaining
    their relative ratios), and all coeffs with priority {p>4} are
    zero. */

void cpk_print_weights(FILE *wr, cpk_wcoeffs_t *wc);
  /* Prints the weight to file {wr}. */

int32_t cpk_weight_cmp(double Wa, double Wb);
  /* Returns 0 if weight {Wa} is equal to {Wb} apart from
    roundoff-like difference; otherwise returns -1 if {Wa < Wb}, +1 if
    {Wa > Wb}. */

#endif
