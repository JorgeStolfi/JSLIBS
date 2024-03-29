#ifndef msm_cand_refine_H
#define msm_cand_refine_H

/* Scale-to-scale mapping and refinement of candidates and candidate lists. */
/* Last edited on 2022-10-20 06:40:53 by stolfi */

#define msm_cand_refine_H_COPYRIGHT \
  "Copyright � 2005  by the State University of Campinas (UNICAMP)" \
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
    int32_t delta,
    int32_t kappa,
    int32_t expand,
    int32_t shrink,
    int32_t maxUnp, 
    msm_rung_step_score_proc_t *step_score,
    bool_t verbose,
    msm_dyn_tableau_t *tb,
    int32_t *n_steps,          
    int32_t *n_entries         
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
    
    The parameters {expand} and {shrink} must be non-negative integers.
    The range of {R}-coordinates of the refined candidate may be longer
    than the original range, by up to {expand} units, or shorter, by up
    to {shrink} units, at each end, if that results in a higher score.
    In particular, if {expand} is zero, the range of {R}-coordinates
    spannned by the candidate will not be allowed to increase. If
    {shrink} is 0, the range of {R}-coordinates will not be allowed to
    decrease; it will be a superset of the original range. If both
    {expand} and {shrink} are zero, therefore, the refined candidate
    will span exactly the same {R}-coordinates as the original.
    
    In any case, the candidate will be confined to the set of valid
    rungs, as defined by {delta} and {kappa}. Thus, the refined
    {R}-range will never exceed the original one by more than {2*(kappa
    + delta)}, in each direction, even if {expand} is larger than that.
    
    The procedure uses the tableau {tb} as a work area, expanding it
    as needed. 
    
    The procedure increments {n_entries} with the number of tableau
    entries that were actually computed and {n_steps} with the number
    of steps that were considered to compute those entries.
    
    If {verbose} is true, the original and refined candidates are
    printed to {stderr}. */

#endif

