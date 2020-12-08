#ifndef msm_cand_vec_H
#define msm_cand_vec_H

/* Tools for lists of candidates. */
/* Last edited on 2017-04-28 10:31:02 by stolfilocal */

#define msm_cand_vec_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_seq_desc.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_dyn.h>

#include <vec.h>

#include <stdint.h>

vec_typedef(msm_cand_vec_t, msm_cand_vec, msm_cand_t);
  /* A {msm_cand_vec_t} is a vector of {msm_cand_t}. */

msm_cand_vec_t msm_cand_vec_get_best_perfect
  ( msm_seq_desc_t *seq0,
    msm_seq_desc_t *seq1,
    msm_rung_step_score_proc_t *step_score,
    int64_t minRungs, 
    double minScore,
    int maxCands
  );
  /* Build a set of initial cands between {*seq0} and {*seq1}, which had better have
    have the same filtering and {estep}. Considers only maximal perfect cands with
    at least {minRungs} rungs and score at least {minScore}. Returns
    only the first {maxCands} cands in order of decreasing score.
    
    The score of a pairing is taken to be be the sum of the
    {step_score} for all its steps, plus the scores of the steps from  
    nowhere to the first rung, and from the last rung to nowhere. */

void msm_cand_vec_throw
  ( int ntr,
    msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    int minlen,
    int maxlen,
    double atomProb,
    double diagProb,
    double skipProb,
    msm_cand_vec_t *cdv,
    int *ncdP
  );
  /* Generates {ntr} random candidates between two abstract sequences
    {seq0,seq1}. The candidates are stored into the vector {*cdv}
    starting at index {*ncdP}, which is incremented with {ntr}.
    
    The number of rungs of each candidte is chosen in the range {minlen..maxlen}
    with uniform probability.
    
    With probability {atomProb}, the pairings will be atomic. 
    
    With probability {diagProb}, the pairings will start and end near
    the diagonal of the rectangle {{0..seq0->size-1}×{0..seq1->size-1}}.
    
    The other parameters are as explained under {msm_cand_throw}. */
 
msm_cand_vec_t msm_cand_vec_map(msm_cand_vec_t *cdv, msm_seq_desc_t *s0_new, msm_seq_desc_t *s1_new);
  /* Applies {msm_cand_map(cd,s0_new,s1_new)} to each element of the candidate 
    vector {cdv}, and returns a vector with the results. 
    
    Note that the pairings of the result may have non-atomic and/or
    non-increasing steps. */

msm_cand_vec_t msm_cand_vec_interpolate(msm_cand_vec_t *cdv);
  /* Applies {msm_cand_interpolate(cd)} to every candidate {cd} in {cdv}. */

msm_cand_vec_t msm_cand_vec_make_increasing(msm_cand_vec_t *cdv, int minIncrEach, int minIncrSum);
  /* Applies {msm_cand_make_increasing(cd,minIncrEach,minIncrSum)}
    to every candidate {cd} in {cdv}. */

msm_cand_vec_t msm_cand_vec_refine
  ( msm_cand_vec_t *cdv,
    msm_seq_desc_t *ap, 
    msm_seq_desc_t *bp,
    int delta,
    int kappa,
    int expand,
    int shrink,
    int maxUnp, 
    msm_rung_step_score_proc_t *step_score, 
    bool_t verbose,
    msm_dyn_tableau_t *tb, 
    int minCover,
    int maxCands,
    double frac
  );
  /* Re-optimizes the candidates {cdv[...]}, which should belong to
    the same two sequences {ap,bp}. The two sequences had better have
    the same filtering and {estep}. The parameters
    {delta,kappa,expand,shrink,maxUnp,step_score,tb}
    are passed to {msm_cand_refine} for each candidate.
    
    Keeps at most {maxCands} candidates, and discards any candidates
    that span less than {minCover} consecutive positions in either
    sequence, or that are (nearly) contained in worse or equivalent
    ones. A candidate {A} is assumed to be contained in another
    candidate {B} if the fraction of rungs of {A} that are in {B} is
    {frac} or more.
    
    The {verbose} flag is passed to each {msm_cand_refine} call. */

void msm_cand_vec_insert
  ( msm_cand_vec_t *cdv,
    int *ncandP,
    int minCover,
    int maxCands,
    msm_cand_t *cd,
    double frac
  );
  /* Assumes that {cdv->e[0..*ncandP-1]} is a set of candidates, sorted
    by decreasing score.  Tries to store into the {cdv} list the pairing {cd}. 
    
    The procedure is a no-op if the pairing is not among the best
    {cd.ne} entries found so far.
    
    The procedure discards {cd} if it spans less than {minCover}
    consecutive positions on either sequence.
    
    If {frac} is positive and no greater than 1, the procedure discards
    any candidate {ca} from {cd} or {cdv} that is mostly contained in
    some other candidate {cb} with better or equal score; where `mostly
    contained' means that the rungs of {ca} that are present in {cb} are
    at least {frac} times the number of rungs of {ca}. */
    
void msm_cand_vec_free(msm_cand_vec_t *cdv);
  /* De-allocates all storage associated with {*cd}; including
    the vector {cdv->e}, but not the record {*cdv} itself. */

/* CANDIDATE LIST I/O */

void msm_cand_vec_write(FILE *wr, msm_cand_vec_t *cdv);
  /* Writes the list of candidates {cdv} to the file {wr}. */

void msm_cand_vec_write_named(msm_cand_vec_t *cdv, char *name, char *tag, char *ext);
  /* Writes the list of pairings {cdv} to the file named "{name}{tag}{ext}". */

msm_cand_vec_t msm_cand_vec_read(FILE *rd);
  /* Reads a list of candidates {cdv} from file {rd}. */

msm_cand_vec_t msm_cand_vec_read_named(char *name, char *tag, char *ext);
  /* Reads a list of candidates {cdv} from file "{name}{tag}{ext}". */

#endif
