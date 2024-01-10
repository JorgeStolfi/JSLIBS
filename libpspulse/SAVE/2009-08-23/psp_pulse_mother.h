/* Univariate polynomial mother pulses */
/* Last edited on 2009-08-23 19:50:48 by stolfi */

#ifndef psp_pulse_mother_H
#define psp_pulse_mother_H

#include <bz_basic.h>
#include <psp_basic.h>
#include <psp_grid.h>
#include <psp_pulse.h>
  
/*
  MOTHER PULSES
  
  Each mother pulse is a bounded-support spline defined on the regular
  infinite unidimensional grid {U} with unit-length intervals.
  
  MOTHER PULSE SUPPORT
  
  The support of a mother pulse {w} always spans a certain positive
  number {msz} of consecutive grid cells (unit intervals) and possibly
  their inferior corner vertices. The number {msz} is the /mother's
  support count/.
  
  By definition, if the support count of a mother pulse is {msz}, the support 
  consists of cells with indices {0..msz-1} and their vertices, whose union 
  is the interval {(0 _ msz)}.
  
  A mother pulse with continuity order {c = -1} is generally undefined
  for integer arguments, but is zero for any argument {x} with {x < 0}
  or {x > msz}. A mother pulse with {c \geq 0} is defined everywhere
  and its value is zero for any argument {x} with {x \leq 0} and
  {x\geq msz}. */

psp_pulse_mother_count_t psp_pulse_mother_count
  ( psp_pulse_kind_t pkind, 
    psp_cont_t c, 
    psp_degree_t g
  );
  /* Number of distinct mother pulses in the family defined by pulse
    kind {pkind}, continuity {c}, and degree {g}. Returns zero iff the
    combination {pkind,c,g} is invalid. */
 
psp_grid_size_t psp_pulse_mother_supp_count
  ( psp_pulse_family_t *fam, 
    psp_pulse_mother_index_t pix
  );
  /* Computes the support count {msz} of the mother spline of 
    family {fam} and mother pulse index {pix}. */

/* MOTHER PULSE EVALUATION

  The following procedures refer to the mother pulse of  
  family {fam}, and mother pulse index {pix}. */

void psp_pulse_mother_eval
  ( psp_pulse_family_t *fam, 
    psp_pulse_mother_index_t pix, 
    double z,
    psp_cont_t ord,
    double *w
  );
  /* Evaluates the mother pulse of family {fam} and mother index {pix},
    and all its derivatives up to order {ord}, at the argument {z}.
    The pulse is specified by its
    
    Returns the result in {w[0..ord]}: the pulse's value in {w[0]},
    and its {k}th derivative in {w[k]}. Beware that {w[k]} may be
    undefined when {k > c} and {z} is an integer. All these values are
    zero if {z} is outside the interval {[0 _ msz]} where 
    {msz = psp_pulse_mother_supp_count(fam,pix)}. */
    
void psp_pulse_mother_to_bezier
  ( psp_pulse_family_t *fam, 
    psp_pulse_mother_index_t pix,
    psp_grid_pos_t ix,
    double *bz
  );
  /* Stores into {bz[0..g]} the Bézier coefficients of the mother
    pulse in the interval {(ix _ ix+1)}. The coefficients will be zero
    if {ix < 0} or {ix >= msz}, where 
    {msz = psp_pulse_mother_supp_count(fam,pix)}. */
  
/*
  Pulses that differ only in grid size and shift are related to
  each other by argument scaling and translation (like wavelets). The
  other four parameters -- kind {pkind}, continuity {c}, degree {g}, and
  mother pulse index {pix} -- define the basic shape of the pulse, and
  thus will be called the /shape parameters/.
  
  All dyadic pulses with the same shape parameters are derived from a
  single /mother pulse/ {w(x)}. The mother pulse is a polynomial
  spline of continuity order {c} and degree {g}, defined on the
  *non-periodic* integer grid -- the regular infinite grid on the real
  line with integer vertices. */

#endif

