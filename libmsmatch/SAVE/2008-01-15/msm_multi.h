#ifndef msm_multi_H
#define msm_multi_H

/* Multiscale DNA matching */
/* Last edited on 2008-01-12 12:13:28 by stolfi */ 

#define msm_multi_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>

typedef void msm_multi_report_proc_t
  ( msm_seq_desc_t *ap, 
    msm_seq_desc_t *bp,
    msm_cand_vec_t *cdv,
    bool_t refined
  );
  /* Type of a procedure that is called by {msm_multi_find_matches} to 
    report the state of the computation at each scale. 
    
    The parameters {ap} and {bp} are the sequences, reduced to the
    same scale (indicated by their {level} field). The parameter {cdv}
    is a list of candidate pairings between the two sequences.
    
    If {refined} is FALSE, {cdv} is the `raw' list of candidates at
    this step (either the initial `brute-force' candidates, or the
    result of mapping and interpolation).
    
    If {refine} is TRUE, {cdv} is the final list of candidates, after
    refinement (if any), sorting, and pruning. */

msm_cand_vec_t msm_multi_find_matches
  ( msm_cand_vec_t *cdvini,
    msm_seq_desc_t a[], 
    msm_seq_desc_t b[],
    int maxMatches,
    int maxLevel,
    int minBases,
    int nwtb0,
    int nwtb1,
    double frac,
    bool_t refine,
    int delta,
    int kappa,
    bool_t shrink,
    int maxunp,
    msm_rung_step_score_proc_t *step_score,
    msm_multi_report_proc_t *report
  );
  /* Finds good pairings between two sequences {a[0]} and {b[0]} by
    multiscale matching up to recursion level {maxLevel}. Assumes that
    {a[r]} and {b[r]} are the sequences {a[0]} and {b[0]} filtered to
    scale {r} (and with the {level} field equal to {r}) for each {r} in
    {0..maxLevel}. Returns the {maxMatches} candidate pairings
    with largest scores, among those that have {minBases} or more bases.
    
    The client must provide the initial candidate list {cdvini}
    for sequences {a[maxLevel]} and {b[maxLevel]}. One can use 
    {msm_cand_get_best_perfect} to obtain such a list. 
    
    If the {refine} flag is true, at each scale {r} each candidate is
    re-optimized and its score is re-computed. Then the candidate list
    is trimmed to have at most {maxMatches*2^r} candidates, after
    discarding those that correspond to less than {minBases} bases in
    the finest scale.
    
    If {refine} is TRUE and {frac >= 0}, the procedure also discards
    any refined candidate {ca} that is mostly contained in some other
    candidate {cb} with better or equal score; where `mostly
    contained' means that the rungs of {ca} that are present in {cb}
    are at least {frac} times the number of rungs of {ca}.
    
    The parameters {nwtb0} and {nwtb1} should be the number of
    elements in the weight tables used for filtering: {nwtb0} from
    level 0 to level 1, {nwtb1} for all the other levels.
    
    The paramters {delta,kappa,shrink,maxunp} are explained under
    {msm_cand_refine}.
    
    If {report} is not NULL, the procedure calls {report(ap, bp,
    cdvmap, cdvfin)} at each level of the matching procedure, from
    {maxLevel-1} down to 0; once for the `raw' candidate list, and
    again for the refined, pruned and sorted list. */

#endif
