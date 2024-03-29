#ifndef msm_cand_H
#define msm_cand_H

/* A candidate is a pairing with a score. */
/* Last edited on 2008-01-12 12:01:09 by stolfi */

#define msm_cand_H_COPYRIGHT \
  "Copyright � 2005  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_seq_desc.h>
#include <msm_pairing.h>

#include <vec.h>

#include <stdint.h>

typedef struct msm_cand_t
  { msm_seq_desc_t seq[2]; /* Descriptors of the two sequences. */
    msm_pairing_t *pr;     /* A partial pairing between the two sequences. */
    double score;          /* Value (goodness) of cand, for sorting. */
  } msm_cand_t;
  /* A candidate pairing between two abstract sequences 
    whose essential attributes are {seq[0],seq[1]}.
    
    Each rung of the pairing {pr} must consist of a valid index {c[0]}
    into {seq[0]}, and a valid index {c[1]} into the {seq[1]}. For
    each {j}, in particular, the first rung must have {c[j]} in the
    range {0..seq[j].nbas-1}. Moreover, if {seq[j].circ} is
    FALSE, this constraint must be satisfied for every rung.
    
    In the current implementation, the pairing must be 1-increasing.
    Note that if {seq[j].circ} is TRUE the pairing may have indices
    {c[j]} that exceed {seq[j].nbas}. Those indices should be taken
    modulo {seq[j].nbas} when indexing into sequence {seq[j]}.
    
    The pairing {pr} itself may be circular. In that case one must
    have {seq[0].circ==seq[1].circ==TRUE}. Moreover, the pairing's
    coordinate period vector {per = msm_pairing_period(pr)} must have
    {per.c[j] % seq[j].nbas == 0} and {per.c[j] > 0} (i.e., it must
    make a positive number of full turns around sequence {seq[j]}).
    
    A circular pairing may cover the same part of {seq[j]} more than
    once. This situation may be useful during multiscale search (to
    represent multiple pre-candidates between two circular sequences),
    and may even make sense for some applications. */

msm_cand_t msm_cand_from_pairing
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp,
    msm_pairing_t *pr,
    double score
  );
  /* Creates a candidate between the two sequences {xp,yp} with
    pairing {pr} and the given {score}. The sequences must have the
    same {level}.
    
    If the pairing {pr} is circular, each coordinate period must be
    positive and a multiple of the corresponding sequence length. The
    two sequences must be circular.
    
    If the pairing {pr} is open, each sequence must be either circular
    or long enough to contain the corresponding endpoint of the last
    rung of {pr}. */

bool_t msm_cand_equivalent(msm_cand_t *cda, msm_cand_t *cdb, bool_t die);
  /* Returns TRUE iff the two candidates have the same sequence
    descriptors and the same rungs. Otherwise, if {die} is true, fails
    with an error message; if {die} is FALSE, returns FALSE silently.
    Does not compare the scores. */

int msm_cand_count_common_rungs(msm_cand_t *ca, msm_cand_t *cb);
  /* Number of rungs that are shared between {ca} and {cb}.
    The two candidates must refer to the same two sequences,
    otherwise the result is 0. Takes circularity of the 
    sequences into account. */

double msm_cand_compute_score
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp,
    msm_pairing_t *pr,
    msm_rung_step_score_proc_t *step_score
  );
  /* Computes a numeric score for the pairing {pr} between the 
    indices of the two sequences {xp} and {yp}.
    
    The score is the sum of {step_score} over all fundamental steps of
    the pairing. These include the loopback step, if the pairing is
    circular; and the none-to-first-rung and the last-rung-to-none
    steps, if it is open. */

/* RANDOM CANDIDATES */

msm_cand_t msm_cand_throw
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp, 
    int len,
    bool_t circp,
    bool_t atomic,
    bool_t diagonal,
    double skipProb
  );
  /* Generates a random candidate between two abstract sequences
    {xp,yp}. The pairing is generated by
    {msm_pairing_throw(xp->nbas,xp->circ,yp->nbas,yp->circ,len,circp,atomic,diagonal,skipProb)}. */

/* CANDIDATE MAPPING */

msm_cand_t msm_cand_map_to_finer(msm_cand_t *cd, int nxnew, int nynew, int nwtb);
  /* Given a candidate {cdold} at some positive scale {r}, returns the
    corresponding candidate {cdnew} at scale {r-1}.
    
    The candidate's pairing is mapped with
    {msm_filter_map_pairing_to_finer}. Assumes that the two sequences
    being paired have {nxnew} and {nynew} elements at the finer scale {r-1},
    and that the filtering used a weight table of {nwtb} entries. The
    candidate's score is preserved. */

msm_cand_t msm_cand_interpolate(msm_cand_t *cd);
  /* Returns a copy of {cd} whose pairing has been interpolated as needed 
    to produce an atomic pairing. */

/* CANDIDATE I/O */

msm_cand_t msm_cand_read(FILE *rd, char *pre, char *suf);
  /* Reads a candidate from {rd}, in the format 
  
      "{pre} {seq0} {seq1} {score} {pairing} {suf}"
  
    The fields {seq1} and {seq2} should be the descriptors of the two
    sequences, in the format accepted by {msm_seq_desc_read} with
    delimiters {"("} and {")"}. 
    
    The {score} should be a float number. 
    
    The {pairing} should be in the format accepted by
    {msm_pairing_read}. 
    
    In general, consecutive fields (and their subfields) must be
    separated by one or more spaces. If the arguments {pre} and {suf}
    are not NULL, they must be present in the input, and must not
    contain embedded blanks or other formatting char. */

void msm_cand_write
  ( FILE *wr, 
    char *pre, 
    msm_cand_t *cd, 
    int idSize, 
    int nameSize, 
    int ixSize, 
    char *suf
  );
  /* Writes the candidate {cd} to the file {wr}, in a format that is
    accepted by {msm_cand_read}.  The {pre} and {suf} strings are
    omitted if NULL.
    
    The numeric IDs and names of sequences are blank-padded to
    {idSize} and {nameSize} characters, respectively. Sequence lengths
    and indices are padded to {ixSize} digits. */

void msm_cand_debug(char *tag, msm_cand_t *cd);
  /* Prints the {tag} and the candidate {cd} to {stderr}. */

/* VALIDATION */

bool_t msm_cand_is_valid(msm_cand_t *cd, bool_t die);
  /* Runs consistency tests on the candidate {cd}. If the tests
    succeed, returns TRUE.  Otherwise, if {die} is true, fails with an
    error message; if {die} is FALSE, returns FALSE silently. */

#endif
