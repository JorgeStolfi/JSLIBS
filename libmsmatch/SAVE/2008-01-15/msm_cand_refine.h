#ifndef msm_cand_refine_H
#define msm_cand_refine_H

/* Scale-to-scale mapping and refinement of candidates and candidate lists. */
/* Last edited on 2008-01-12 12:12:59 by stolfi */

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
    bool_t shrink,
    int maxunp, 
    msm_rung_step_score_proc_t *step_score,
    msm_dyn_tableau_t *tb 
  );
  /* Finds a pairing with maximum score between two sequences {ap} and
    {bp} within the band of half-width {delta} around the candidate
    {cdold}, with ends possibly extended by {kappa} along either or
    both sequences, in both directions.
    
    The candidate {cdold} should contain a strictly increasing pairing
    between the two sequences {ap} and {bp}. That pairing is
    implicitly interpolated as necessary to make it atomic and
    strictly increasing.
    
    The new pairing is constrained to use only /valid/ rungs and /valid/
    steps.
    
    The set of valid rungs is defined as the subset of the integer XY
    grid which is the union of squares of side {2*delta+1} centered on
    rungs of the old pairing, united with squares of side {2*kappa+1}
    whose maximum (resp. minimum) corner [check !!!] are the first
    (resp. last) rung of that pairing.
    
    A step {f-->g} with X and Y displacements {dx,dy} is considered
    valid if and only if (0) both rungs are valid; (1) the step is
    strictly increasing (that is, {dx > 1} and {dy > 1}); and (2) the
    step leaves at most {maxunp} elements unpaired, that is, 
    {dx+dy-2 <= maxunp}.
    
    The current implementation supports only open pairings.
    
    If {shrink} is TRUE, the procedure allows the candidate to shrink
    arbitrarily. If {shrink} is false, the end-rungs of the result are
    constrained to be inside the squares of side {2*kappa+1} centered
    on the original end-rungs. Uses the tableau {tb} as a work area,
    expanding it as needed. */

#endif

