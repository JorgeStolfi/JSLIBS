#ifndef msm_cand_H
#define msm_cand_H

/* A candidate is two sequence descriptors, a pairing between them, and a score. */
/* Last edited on 2022-10-20 06:41:16 by stolfi */

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
  /* An {msm_cand_t} describes a candidate matching between 
    elements of two abstract sequences whose essential attributes are
    {seq[0],seq[1]}.
    
    The pairing {pr} may be empty. If not empty, its rungs must
    be pairs {(i0,i1)} with {i0} in {0..seq[0].size-1} and {i1} in
    {0..seq[1].size-1}. Also, every step of {pr} must be non-decreasing
    on both sides and non-stuttering (must be strictly increasing on at
    least one side). Note that not every valid pairing is valid for a
    candidate.  Some applications may impose additional constraints
    on the pairing. */

msm_cand_t msm_cand_from_pairing
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    msm_pairing_t *pr,
    double score
  );
  /* Creates a candidate between the two sequences {seq0,seq1} with
    pairing {pr} and the given {score}.  
    
    The pairing {pr} may be empty. If not empty, its rungs of {pr}
    must be pairs in {{0..seq0.size-1}�{0..seq1.size-1}}, and its steps
    must be non-decreasing on both sides and non-stuttering
    (must be strictly increasing on at least one side). */

void msm_cand_free(msm_cand_t *cd);
  /* De-allocates all storage associated with {*cd}, except the
    record {*cd} itself. */

bool_t msm_cand_equivalent(msm_cand_t *cda, msm_cand_t *cdb, bool_t die);
  /* Returns TRUE iff the two candidates have the same sequence
    descriptors and the same rungs. Otherwise, if {die} is true, fails
    with an error message; if {die} is FALSE, returns FALSE silently.
    Does not compare the scores. */

int32_t msm_cand_count_common_rungs(msm_cand_t *ca, msm_cand_t *cb);
  /* Number of rungs that are shared between {ca} and {cb}.
    The two candidates must refer to the same two sequences,
    otherwise the result is 0. */

double msm_cand_compute_score
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    msm_pairing_t *pr,
    msm_rung_step_score_proc_t *step_score,
    bool_t ini_step,
    bool_t fin_step
  );
  /* Computes a numeric score for the pairing {pr} between the 
    indices of the two sequences {seq0} and {seq1} (which had better
    have the same filtering and {estep}).
    
    The score is the sum of {step_score} over all steps of the pairing.
    If {ini_step} is true, adds the score of the nowhere-to-first-rung.
    If [fin_step} is true, adds the score of the last-rung-to-nowhere
    steps. */

/* RANDOM CANDIDATES */

msm_cand_t msm_cand_throw
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1, 
    int32_t len,
    bool_t atomic,
    bool_t diagonal,
    double skipProb
  );
  /* Generates a random candidate between two abstract sequences
    {seq0,seq1}. The pairing is generated by {msm_pairing_throw(seq0->size,
    seq1->size, len,atomic,diagonal,skipProb)}. */

msm_cand_t msm_cand_sub_throw(msm_cand_t *cd, int32_t len);
  /* Extracts from candidate {cd} a candidate between the same
    sequences whose pairing is a set of consecutive rungs
    of {cd}'s pairing, as generated by
    {msm_pairing_sub_throw(cd->pr,len)}. */

/* CANDIDATE MAPPING */

msm_cand_t msm_cand_map(msm_cand_t *cd,  msm_seq_desc_t *s0_new, msm_seq_desc_t *s1_new);
  /* Maps every rung of the candidate {cd} so that it relates the
    corresponding indices of {s0_new} and {s1_new}; see
    {msm_rung_vec_map_coords} for details.
    
    The candidate's score is just copied, not recomputed. Note that the
    pairing of the result may have non-atomic steps and/or steps 
    that do not increase on one side, or even on both sides. */

msm_cand_t msm_cand_interpolate(msm_cand_t *cd);
  /* Returns a copy of {cd} whose pairing {p} has been 
    processed with {msm_pairing_interpolate(p)}. */

msm_cand_t msm_cand_make_increasing(msm_cand_t *cd, int32_t minIncrEach, int32_t minIncrSum);
  /* Returns a copy of {cd} whose pairing {p} has been processed with
    {msm_pairing_make_increasing(p,minIncrEach,minIncrSum)}. */

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
    contain embedded blanks or other formatting char.
    
    Note that the format of the sequence descriptor changed on 2013-10-16. */

void msm_cand_write
  ( FILE *wr, 
    char *pre, 
    msm_cand_t *cd, 
    int32_t idSize, 
    int32_t nameSize, 
    int32_t ixSize, 
    char *suf,
    bool_t writePairing
  );
  /* Writes the candidate {cd} to the file {wr}, in a format that is
    accepted by {msm_cand_read}.  The {pre} and {suf} strings are
    omitted if NULL.
    
    The numeric IDs and names of sequences are blank-padded to
    {idSize} and {nameSize} characters, respectively. Sequence lengths
    and indices are padded to {ixSize} digits.
    
    if {writePairing} is false, omits the pairing. */

void msm_cand_debug(char *tag, msm_cand_t *cd);
  /* Prints the {tag} and the candidate {cd} to {stderr}. */

/* VALIDATION */

bool_t msm_cand_pairing_is_valid(msm_pairing_t *pr, int32_t n0, int32_t n1, int32_t minIncrEach, bool_t atomic, bool_t die);
  /* Checks whether {pr} is an acceptable pairing for a candidate.
    match between sequences with {n0} and {n1} samples.
    
    Specifically, {pr} must satisfy {msm_pairing_is_valid(pr)},
    {msm_pairing_is_increasing(pr,imax(0,minIncrEach),1), and (if
    {atomic} is true) {msm_pairing_is_atomic(p)}.  Also,
    every rung must be a pair {(i0,i1)} with {i0} in {0..n0-1}
    and {i1} in {0..n1-1}.
    
    If the tests succeed, returns TRUE. Otherwise, if {die} is true,
    fails with an error message; if {die} is FALSE, returns FALSE
    silently. */

bool_t msm_cand_is_valid(msm_cand_t *cd, int32_t minIncrEach, bool_t atomic, bool_t die);
  /* Runs consistency tests on the candidate {cd}.
  
    In particular, requires the pairing {cd.pr} to pass
    {msm_cand_pairing_is_valid(pr,n0,n1,minIncrEach,atomic)}, where
    {n0,n1} are the sizes of the sequences {cd.seq[0],cd.seq[1]}.
    
    If the tests succeed, returns TRUE. Otherwise, if {die} is true,
    fails with an error message; if {die} is FALSE, returns FALSE
    silently. */

#endif
