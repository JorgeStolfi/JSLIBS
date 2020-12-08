#ifndef msm_seq_desc_H
#define msm_seq_desc_H

/* Abstract sequence descriptors. */
/* Last edited on 2017-04-28 16:56:15 by stolfilocal */

#define msm_seq_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <stdint.h>

#include <vec.h>

#include <msm_basic.h>

/* SEQUENCE DESCRIPTOR 

  The main data type {msm_seq_desc_t} below is a descriptor for a finite
  sequence of /sample values/ over some /sample space/ {V} that depends
  on the application.
  
  The sequence is assumed to be the result of upsampling or downsampling
  some /original sequence/ (e.g. RNA/DNA, digitized sound, digitized
  fragment outline) by some power of two. 
  
  The descriptor specifies the {size} (number of samples) in the
  sequence. The valid /sample indices/ are {0..size-1}.
  
  This module does not provide for storage of the sample
  values, which may even be generated on-the-fly by an analytic
  formula or interpolation of some coarser sequence.
  
  The descriptor also specifies a correspondence between the valid
  sample indices and the sample indices of the original sequence,
  assuming that it is a simple affine map. (If the correspondence is not
  an affine map, it must be represented in some other way). */

typedef int msm_seq_id_t; 
  /* A {msm_seq_id_t} is an internal identifier for an original
    data sequence, e.g. an index into a table of biosequences. */

#define msm_seq_id_none (-1)
  /* A {msm_seq_id_t} value that means "no original sequence". */
  
#define msm_seq_desc_estep_MAX 24
  /* Should not need more than this much upsampling/downsampling.
    This limit must be less than 30 so that {2^{abs(estep)}}
    can be safely computed in 32 bits. */

typedef struct msm_seq_desc_t 
  { /* Original sequence identification: */
    msm_seq_id_t id;  /* Internal identifier of the sequence. */
    char *name;       /* External identifier of the sequence. */
    bool_t rev;       /* TRUE if the original sequence was reversed. */
    /* Sequence size: */
    int32_t size;     /* Number of valid samples in the sequence. */
    /* Relation to the original sequence: */
    int8_t estep;     /* Each sampling step is {2^estep} steps of original seq; may be negative. */
    int32_t skip;     /* First sample matches sample {skip*(2^estep)} of original seq. */
  } msm_seq_desc_t;
  /* The {id} and/or the {name} fields identify the original
    sequence from which this sequence was derived.  
    
    The {rev} flag, if true, indicates that the sequence is derived from
    the reversal of the original sequence indicated by {id} and {name}.
    (For some applications, e.g. DNA/RNA or fragment outline matchings, this flag
    also indicates complementation of sample values.)
    
    The valid sample indices for this sequence are {0..size-1}.
    
    If {size} is zero,{skip} and {estep} are irrelevant.
    
    If {size} is nonzero, then sample index {i} of this sequence
    corresponds to sample {j = (i + skip)*{2^estep}} of the original
    sequence (or of the reversed sequence, if {rev} is false).
    
    Note that some integer indices of {s} may correspond to fractional
    indices of the oriiginal sequence, if {estep} is negative.
    
    The absolute value of {estep} must not exceed {msm_seq_desc_estep_MAX}. */

msm_seq_desc_t msm_seq_desc_make
  ( msm_seq_id_t id,  /* Internal identifier of the sequence. */
    char *name,       /* External identifier of the sequence. */
    bool_t rev,       /* TRUE if the original sequence was reversed. */
    int32_t size,     /* Number of distinct values in the sequence. */
    int8_t estep,     /* Each sampling step is {2^estep} samples of original; may be negative. */
    int32_t skip      /* First sample matches sample {skip*(2^estep)} of original seq. */
  );
  /* Assembles an {msm_seq_desc_t} from the given fields.
    Fails if {size} is negative, or the parameters are not
    valid as specified under {msm_seq_desc_t}. */

/* DESCRIPTOR COMPARISON

  Each procedures in this section compares two descriptors {sa,sb} in
  various ways. If the comparison succeeds, the precedure returns TRUE.
  Otherwise, if {die} is true, the procedure fails with an error
  message; if {die} is FALSE, it returns FALSE silently. The {name}
  strings are compared with {strcmp}. */

bool_t msm_seq_desc_same_orig_seq(msm_seq_desc_t *sa, msm_seq_desc_t *sb, bool_t die);
  /* Returns TRUE if the sequence descriptors {sa} and {sb} are derived
    from the same original sequence (that is, if they have the
    same {name}, {id}, and {rev} fields). */

bool_t msm_seq_desc_equal(msm_seq_desc_t *sa, msm_seq_desc_t *sb, bool_t die);
  /* Returns TRUE if the sequence descriptors {sa} and {sb} are derived
    from the same original sequence and have the same sample positions
    relative to it (that is, if they have exactly the same fields.
    
    The fields {estep} and {skip} are compared only if the {size} fields
    match and are positive. */

/* DESCRIPTOR MANIPULATION */
 
msm_seq_desc_t msm_seq_desc_trim(msm_seq_desc_t *s, int32_t itrim, int32_t ftrim);
  /* Returns a descriptor for the sequence {t} that would be
    obtained from sequence {s} by discarding the first {itrim}
    samples and the last {ftrim} samples.
    
    The arguments {itrim} and {ftrim} may be zero, or negative to
    account for extrapolation of {s}.
    
    The procedure copies the fields {t.id}, {t.name}, {t.rev}, and
    {t.estep} from {s}, and computes the proper values for {t.size} and
    {t.skip}. If the sequence {s} has fewer than {itrim + ftrim}
    samples, then {t.size} and {t.skip} are set to zero. */

msm_seq_desc_t msm_seq_desc_resample(msm_seq_desc_t *s, int8_t ek);
  /* Returns a descriptor for the sequence {t} that would be obtained by
    resampling a sequence with descriptor {s} with one new sample
    every {2^ek} samples of {s}.
    
    In all cases, the result {t} has the same fields {id}, {name}, and {rev}
    as {s}, and has {t.estep = s.estep + ek}.
    
    If {ek} is zero, the result {t} is the same as {s}.
    
    If {s} has zero samples, the result {t} too has zero samples,
    and {t.skip} is set to zero.
    
    Otherwise, the {skip} and {size} of {t} are computed so that the
    first and last samples of {t} are aligned with samples of {s}.
    
    Namely, if {ek} is positive, the procedure assumes that {t} is
    obtained by taking one every {m} samples of {s}. The sequence {s} is
    first implicitly trimmed at both ends, if needed, so that {s.skip}
    and {s.size-1} are multiples of {m}.
    
    If {ek} is negative, the procedure assumes that {t} is obtained
    by inserting additional {m-1} samples between every two consecutive
    samples of {s}, where {m = 2^{-ek}}. In this case, the first and
    last sample positions of {t} will correspond to the first and last
    ones of {s}. That is, {t.size} will be {(s.size - 1)*m + 1}, and
    {t.skip} will be {s.skip*m}. */ 

msm_seq_desc_t msm_seq_desc_filter(msm_seq_desc_t *sd, int nw, int8_t ek);
  /* Returns an abstract sequence descriptor for the result of
    filtering an abstract sequence with descriptor
    {sd}, with a filter whose table has {nw} entries, and resampling the result
    with step {2^ek} times the old step.
    
    The output descriptor {fd} will have the same attributes {id,name,rev} as
    {sd}, and its {estep} will be {sd.estep+ek}.  The {size}
    and {skip} will be adjusted so that the first and last samples 
    of {fd} can be computed by convolution with that kernel fitting entirely 
    inside {sd} and centered over a sample of the original biosequence. */

/* INDEX MAPPING BETWEEN SUCCESSIVE SCALES OF RESOLUTION */

double msm_seq_desc_map_index_to_orig_seq(double is, msm_seq_desc_t *s);
  /* Maps a sample index {ix} in sequence {s} to the
    corresponding sample index {ix0} of the original  
    sequence {s0}.  
    
    The given index {ix} is not required to be in {0..s.size-1}
    but must be an integer or a dyadic fraction.
    
    The result {ix0} may be outside {0..s0.size-1}, and may be a dyadic
    fraction if {ix} is a fraction or {s.estep} is negative. Even
    then, the the result will be exact. */

double msm_seq_desc_map_index_from_orig_seq(double ix0, msm_seq_desc_t *s);
  /* Maps a sample index {ix0} in some original sequence {s0}
    to the corresponding sample index in sequence {s}, which
    is assumed to be a filtered and resampled version of {s0}.
    
    The given index {ix0} is not required to be in {0..s0.size-1}
    and must be an integer or a dyadic fraction.
    
    The result {ix} may be outside {0..s.size-1} and may be a dyadic
    fraction, if {ix0} is a fraction or {s.estep} is positive. Even
    then, the the result will be exact. */

double msm_seq_desc_map_index(double ixa, msm_seq_desc_t *sa, msm_seq_desc_t *sb);
  /* Maps a sample index {ixa} in a sequence {sa} to the
    corresponding index {ixb} in another version {sb} of the same 
    original sequence.  
    
    The two sequences must have the same {id}, {name} and {rev} fields.  
    The given index {ixa} must be an integer or a dyadic fraction.
    The result {ixb} will be  an integer or a dyadic fraction.
    This operation is exact in any case. */

/* RANDOM SEGMENTS */

void msm_seq_desc_throw_segment(int32_t size, int32_t n, int32_t *iniP, int32_t *finP);
  /* Selects initial and final sample indices {*ini,*fin} for a random
    segment spanning {n} positions, taken from a sequence with {size}
    positions. */ 

/* PRINTOUT */

void msm_seq_desc_write
  ( FILE *wr, 
    char *pre, 
    msm_seq_desc_t *s, 
    int idSize, 
    int nameSize, 
    int indexSize,
    char *suf
  );
  /* Writes to {wr} the sequence descriptor {s}, in a format
    compatible with {msm_seq_desc_read}; namely,
    
     "{pre} {s.id} {s.name} {s.rev} {s.size} {s.estep} {s.skip} {suf}". 
      
    The delimiters {pre} and {suf} will be omitted if NULL. The fields
    {s.id} and {s.name} will be blank-padded to {idSize} and
    {nameSize} chars.  The field {s.rev} is "F" for false or "T"
    for true.  The field {s.estep} is zero or has a sign "+" or "-".
    
    The fields {s.size} and {s.skip} will be blank-padded to {indexSize} chars.
    
    The format changed on 2013-10-16. Formerly it was
      
     "{pre} {s.id} {s.name} {s.level} {s.size} {s.circ} {suf}"
     
    where {s.level} was a non-negative `filtering level'
    (similar but not equivalent to {s.estep}) and {s.circ} 
    was a circularity flag. */

msm_seq_desc_t msm_seq_desc_read(FILE *rd, char *pre, char *suf);
  /* Parses from {rd} a sequence descriptor, in the format 
  
      "{pre} {id} {name} {rev} {size} {estep} {skip} {suf}". 
      
    If the arguments {pre} and {suf} are not NULL, they must be
    present in the input. In general, consecutive fields should be
    separated by one or more spaces. The {name} field and the
    {pre,suf} arguments must not contain embedded spaces or other
    formatting chars. */

bool_t msm_seq_desc_is_valid(msm_seq_desc_t *s, bool_t die);
  /* Performs several consistency checks in {s}. If it is valid (even if
    empty), returns TRUE. Otherwise, fails with an error message if
    {die} is true, returns FALSE silently if {die} is FALSE. */

#endif
