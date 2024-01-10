#ifndef msm_cand_vec_H
#define msm_cand_vec_H

/* Tools for lists of candidates. */
/* Last edited on 2008-01-12 12:13:08 by stolfi */

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
  ( msm_seq_desc_t *xp,
    msm_seq_desc_t *yp,
    msm_rung_step_score_proc_t *step_score,
    int minBases, 
    double minScore,
    int maxCands,
    double frac
  );
  /* Build a set of initial cands between {*xp} and {*yp}, which must
    have the same {level}. Considers only maximal perfect cands with
    at least {minBases} bases and score at least {minScore}. Returns
    only the best {maxCands} cands in order of increasing score,
    excluding candidates that are contained in better ones at the
    {frac} level of containment.
    
    The score of a pairing is taken to be be the sum of the
    {step_score} for all its steps, including loopback step if the
    pairing is circular, and none-to-first-rung and last-rung-to-none
    if it is open. */


msm_cand_vec_t msm_cand_vec_throw
  ( int ncd,
    msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp,
    int minlen,
    int maxlen,
    double circProb,
    double atomProb,
    double diagProb,
    double skipProb
  );
  /* Generates a vector with {ncd} random candidates between two
    abstract sequences {xp,yp}, which must have the same {level}. 
    
    The candidate lengths are chosen in the range {minlen..maxlen}
    with uniform probability.
    
    If both sequences are circular, the candidates will have circular
    pairings with probebility {circProb}; in that case {maxlen} is
    ignored. Otherwise all candidates will have open pairings.
    
    With probability {atomProb}, the pairings will be 1-atomic. 
    
    With probability {diagProb}, the pairings will start and end near
    the diagonal of the rectangle {{0..xp->nbas-1}×{0..yp->nbas-1}}.
    
    The other parameters are as explained under {msm_cand_throw}. */
 
msm_cand_vec_t msm_cand_vec_map_to_finer(msm_cand_vec_t *cdv, int nxnew, int nynew, int nwtb);
  /* Applies {msm_cand_map_to_finer(cd,nxnew,nynew,nwtb)} to
    each element of the candidate vector {cdv}, and returns a vector
    with the results. */

msm_cand_vec_t msm_cand_vec_interpolate(msm_cand_vec_t *cdv);
  /* Applies {msm_cand_interpolate} to all candidates in {cdv}. */

msm_cand_vec_t msm_cand_vec_refine
  ( msm_cand_vec_t *cdv,
    msm_seq_desc_t *ap, 
    msm_seq_desc_t *bp,
    int delta,
    int kappa,
    bool_t shrink,
    int maxunp, 
    msm_rung_step_score_proc_t *step_score, 
    msm_dyn_tableau_t *tb, 
    int maxCands,
    double frac
  );
  /* Re-optimizes the candidates {cdv[...]}, which should belong to
    the same two sequences {ap,bp}. The two sequences must have the
    same {leve}. The parameters{delta,kappa,shrink,maxunp,step_score,tb}
    are passed to {msm_cand_refine} for each candidate.
    
    Keeps at most {maxCands} candidates, and discards any candidates 
    that are (nearly) contained in worse or equivalent ones.  The 
    parameter {frac} specifies the fraction of rung containment
    that means candidate containment. */

void msm_cand_vec_insert
  ( msm_cand_vec_t *cdv,
    int *ncandP,
    int maxCands,
    msm_cand_t *cd,
    double frac
  );
  /* Assumes that {cdv->el[0..*ncandP-1]} is a set of candidates, sorted
    by decreasing score.  Tries to store into the {cdv} list the pairing {cd}. 
    
    The procedure is a no-op if the pairing is not among the best
    {cd.nel} entries found so far.
    
    If {frac >= 0} the procedure discards any candidate {ca} from {cd}
    or {cdv} that is mostly contained in some other candidate {cb}
    with better or equal score; where `mostly contained' means that
    the rungs of {ca} that are present in {cb} are at least {frac}
    times the number of rungs of {ca}. */
    
/* CANDIDATE LIST I/O */

void msm_cand_vec_write(FILE *wr, msm_cand_vec_t *cdv);
  /* Writes the list of candidates {cdv} to the file {wr}. */

void msm_cand_vec_write_named(msm_cand_vec_t *cdv, char *name, char *tag);
  /* Writes the list of pairings {cdv} to the file named "{name}{tag}.cdv". */

msm_cand_vec_t msm_cand_vec_read(FILE *rd);
  /* Reads a list of candidates {cdv} from file {rd}. */

msm_cand_vec_t msm_cand_vec_read_named(char *name, char *tag);
  /* Reads a list of candidates {cdv} from file "{name}{tag}.cdv". */

#endif
