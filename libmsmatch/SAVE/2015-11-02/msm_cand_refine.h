#ifndef msm_cand_refine_H
#define msm_cand_refine_H

/* Scale-to-scale mapping and refinement of candidates and candidate lists. */
/* Last edited on 2013-10-22 06:58:22 by stolfilocal */

#define msm_cand_refine_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_pairing.h>
#include <msm_dyn.h>
#include <msm_cand.h>

#include <stdint.h>

/* CANDIDATE REFINEMET */

msm_cand_t msm_cand_refine
  ( msm_cand_t *cdold, 
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp, 
    msm_rung_step_score_proc_t *step_score,
    bool_t verbose,
    msm_dyn_tableau_t *tb,
    int *n_steps,          
    int *n_entries         
  );
  /* Finds a pairing with maximum score between two sequences {s0} and
    {s1} within the band of half-width {delta} around the candidate
    {cdold}, with ends possibly extended by {kappa} along either or
    both sequences, in both directions.
    
    The candidate {cdold} should contain a strictly increasing pairing
    between the two sequences {s0} and {s1}. [!!! Check whether this
    description matches the code: !!!] That pairing is implicitly
    interpolated as necessary to make it atomic and strictly increasing.
    
    The new pairing is constrained to use only /valid/ rungs and /valid/
    steps. The set of valid rungs is defined as the subset of the integer XY
    grid which is the union of squares of side {2*delta+1} centered on
    rungs of the old pairing, united with squares of side {2*kappa+1}
    whose maximum (resp. minimum) corner [check !!!] are the first
    (resp. last) rung of that pairing.
    
    A step {f-->g} with displacements {d0,d1} on sides 0 and 1 is considered
    valid if and only if (0) both rungs are valid; (1) the step is
    strictly increasing (that is, {d0 > 1} and {d1 > 1}); and (2) the
    step leaves at most {maxUnp} elements unpaired on either side, that is, 
    {d0 = 1 & d1 <= 1 +  maxUnp} or {d1 = 1 & d0 <= 1 +  maxUnp} .
    
    If {expand} is FALSE, the range of {R}-coordinates spannned by the
    candidate will not be allowed to increase; it will be a subset of
    the original range. If {expand} is TRUE, the {R}-range may expand by
    up to {2*(kappa + delta)} in either direction.
    
    If {shrink} is FALSE, the range of {R}-coordinates will not be
    allowed to decrease; it will be a superset of the original range.
    If {shrink} is TRUE, the refined candidate may be arbitrarily
    shorter than the original, depending on the scores.
    
    If both {expand} and {shrink} are TRUE, therefore, the refined candidate
    will span exactly the same {R}-coordinates as the original. In any case,
    the candidate will be confined to the set of valid rungs, as defined
    by {delta} and {kappa}. 
    
    The procedure uses the tableau {tb} as a work area, expanding it
    as needed. 
    
    The procedure increments {n_entries} with the number of tableau
    entries that were actually computed and {n_steps} with the number
    of steps that were considered to compute those entries.
    
    If {verbose} is true, the original and refined candidates are
    printed to {stderr}. */

#endif

