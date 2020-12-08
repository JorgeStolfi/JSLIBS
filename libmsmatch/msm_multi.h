#ifndef msm_multi_H
#define msm_multi_H

/* Multiscale DNA matching */
/* Last edited on 2017-04-28 10:31:14 by stolfilocal */ 

#define msm_multi_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>

typedef void msm_multi_report_proc_t
  ( int level,
    msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    msm_cand_vec_t *cdv,
    bool_t refined
  );
  /* Type of a procedure that is called by {msm_multi_find_matches} to 
    report the state of the computation at each scale. 
    
    The parameters {seq0} and {seq1} are the sequences, from the 
    same filtering {level}. The parameter {cdv}
    is a list of candidate pairings between the two sequences.
    
    If {refined} is FALSE, {cdv} is the `raw' list of candidates at
    this step (either the initial `brute-force' candidates, or the
    result of mapping and interpolation).
    
    If {refine} is TRUE, {cdv} is the final list of candidates, after
    refinement (if any), sorting, and pruning. */

msm_cand_vec_t msm_multi_find_matches
  ( msm_cand_vec_t *cdvini,
    msm_seq_desc_t sv0[], 
    msm_seq_desc_t sv1[],
    int maxMatches,
    int maxLevel,
    int minSamples[],
    double frac,
    bool_t refine,
    int delta,
    int kappa,
    int expand,
    int shrink,
    int maxUnp,
    msm_rung_step_score_proc_t *step_score,
    msm_multi_report_proc_t *report,
    bool_t verbose
  );
  /* Finds good pairings between two sequences {sv0[0]} and {sv1[0]} by
    multiscale matching up to recursion level {maxLevel}. 
    
    The procedure assumes that {sv0[r]} and {sv1[r]} are the sequences
    {sv0[0]} and {sv1[0]} filtered to scale {r} for each {r} in
    {0..maxLevel}. It does not assume any relationship between {r} ad
    the {.estep} field of {sv0[r],sv1[r]}.

    The procedure returns the {maxMatches} candidate pairings with
    largest scores, among those that span at least {minSamples[0]}
    positions in both sequences.
    
    The client must provide the initial candidate list {cdvini}
    for sequences {sv0[maxLevel]} and {sv1[maxLevel]}. One can use 
    {msm_cand_get_best_perfect} to obtain such a list. 
    
    If the {refine} flag is true, at each scale {r} each candidate is
    re-optimized and its score is re-computed. Then the candidate list
    is trimmed to have at most {maxMatches*2^r} candidates, after
    discarding those that cover less than {minSamples[level]} positions
    on either sequence at that level.
    
    If {refine} is TRUE and {frac >= 0}, the procedure also discards
    any refined candidate {ca} that is mostly contained in some other
    candidate {cb} with better or equal score; where `mostly
    contained' means that the rungs of {ca} that are present in {cb}
    are at least {frac} times the number of rungs of {ca}.
    
    The parameters {delta,kappa,expand,shrink,maxUnp} are explained under
    {msm_cand_refine}.
    
    If {report} is not NULL, the procedure calls {report(seq0, seq1,
    cdvmap, cdvfin)} at each level of the matching procedure, from
    {maxLevel-1} down to 0; once for the `raw' candidate list, and
    again for the refined, pruned and sorted list.
    
    The {verbose} flag requests printing to {stderr} of the original
    and refined candidates at each level. */

#endif
