#ifndef msm_pairing_H
#define msm_pairing_H

/* Pairings between elements of two index ranges */
/* Last edited on 2022-10-20 06:37:43 by stolfi */

#define msm_pairing_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#include <msm_basic.h>
#include <msm_rung.h>

#include <sign.h>
#include <vec.h>

#include <stdint.h>

/* PAIRINGS

  A /pairing/ {pr} is a sequence of one or more rungs {pr[i]}, for some set
  of integer /rung indices/ {i}. The element {pr[i]} is defined only for rung
  indices {i} in {0..ng-1}, for some integer {ng >= 0}.
  
  Every step in a pairing must be non-decreasing on both sides.
  A pairing however may repeat consecutive rungs. */
  
typedef struct msm_pairing_t msm_pairing_t;

int32_t msm_pairing_num_rungs(msm_pairing_t *pr);
  /*Returns the number {ng} of rungs in {pr} */

int32_t msm_pairing_span(msm_pairing_t *pr, int32_t j);
  /* Returns the number {np} of consecutive positions in sequence {j}
    that are spanned by the rungs of pairing {pr}. The result
    is the difference between the {j} components of the first and last
    rungs, plus one. If {j} is -1 returns the sum of both spans.  */

int32_t msm_pairing_sub_span(msm_pairing_t *pr, int32_t ini, int32_t fin, int32_t j);
  /* Returns the number {np} of consecutive positions in sequence {j}
    that are spanned by the rungs of pairing {pr}, between rung number
    {ini} and rung number {fin} inclusive. If {j} is negative, returns
    the sum of the spans on the two sequences.
    
    If {ini > fin}, returns 0. Otherwise, the pairing must be non-empty and
    {ini} and {fin} must be in the range {0..ng-1} where
    {ng=msm_pairing_num_rungs(pr)}. */

msm_rung_t msm_pairing_get_rung(msm_pairing_t *pr, int32_t i);
  /* Obtains the rung with index {i} in the pairing {*pr}. 
    The index {i} must be in the range {0..ng-1}, where {ng} is
    the number of rungs. */

bool_t msm_pairing_equal(msm_pairing_t *pa, msm_pairing_t *pb, bool_t die);
  /* Returns TRUE if the pairings {pa} and {pb} are equal, even with the same 
    repeated rungs.  Otherwise, if {die} is true, fails with an error message; if {die}
    is FALSE, returns FALSE silently. */
  
msm_rung_vec_t msm_pairing_get_rungs(msm_pairing_t *pr);
  /* Returns a new rung vector with a copy of all the rungs of {pr}. */

msm_pairing_t *msm_pairing_empty(void);
  /* Creates an empty pairing. */

msm_pairing_t *msm_pairing_copy(msm_pairing_t *pr);
  /* Creates a storage-independent copy of {*pr}. */
    
msm_pairing_t *msm_pairing_sub_copy(msm_pairing_t *pr, int32_t ini, int32_t fin);
  /* Creates a pairing {s} that is storage-independent from {*pr} and
    contains the rungs between rung number {ini}, inclusive, and rung
    number {fin}, inclusive. 
    If {ini > fin} the result is an empty pairing. Otherwise,
    the procedure requires {0 <= ini <= fin < ng-1},
    where {ng=msm_pairing_num_rungs(pr)}, and the result is non-empty. */
    
msm_pairing_t *msm_pairing_trim_to_span(msm_pairing_t *pr, int32_t j, int32_t span, int32_t i, sign_t dir);
  /* Creates a copy of a segment of {pr} that covers {span} positions on sequence {j}
    (either 0 or 1 for one side of the pairing, or -1 for the sum of both sides).
    The pairing {pr} must be non-empty and {span} must be positive.

    If {dir} is {+1}, the segment will start at rung number {i} of {pr}.
    If {dir} is {-1}, the segment will end at that rung. If {dir} is
    {00}, the segment will include that rung and extend evenly on both
    sides.

    In any case, the returned segment {ps} will be a maximal subset of
    consecutive rungs of {pr} that satisfies {msm_pairing_span(ps,j) <=
    span} and is contained in the original pairing. */

msm_pairing_t *msm_pairing_from_rung_vec(msm_rung_vec_t *gv);
  /* Creates a pairing from a rung vector {gv}.
    
    The pairing will ave {ng = gv.ne} rungs, namely {gv.e[0..ng-1]}. The
    steps must be non-decreasing on both sides, but may be stuttering.
    The order and multiplicity of the rungs is preserved. NOTE: The
    pairing will share the vector {gv->e}. */

void msm_pairing_free_rungs(msm_pairing_t *pr);
  /* Reclaims the internal storage (rung vector) used by {pr}, if any.
    Does not free {*pr} itself. */

void msm_pairing_free(msm_pairing_t *pr);
  /* Reclaims the pairings {*pr}, including the descriptor {*pr} and any
    internal storage (rung vector). */

void msm_pairing_multi_free(msm_pairing_t *pf[], int32_t maxLevel);
  /* Reclaims all the pairings {*(pf[0..maxLevel])}. */

/* PERFECT PAIRINGS
  
   A pairing is /perfect/ if it is non-empty and 
   the increment between any two consecutive rungs is {(+1,+1)}. */

bool_t msm_pairing_is_perfect(msm_pairing_t *pr, bool_t die);
  /* TRUE iff {pr} is a perfect pairing. Otherwise, if {die} is true,
    fails with an error message, else returns FALSE silently. */

msm_pairing_t *msm_pairing_perfect(msm_rung_t gini, int32_t ng);
  /* Creates a perfect pairing with rungs {gini + (k,k)} for {k} in {0..ng-1}. */

msm_rung_vec_t msm_pairing_collect_rungs(msm_pairing_t *pr, int32_t n0, int32_t n1);
  /* Returns the set of all distinct rungs of {pr}, namely rungs  
    {pr[0..ng-1]} where {ng = msm_pairing_num_rungs(pr)}.
    
    The parameters {n0} and {n1} must be non-negative. If {n0} is
    positive, the 0-side index of every rung are reduced modulo {n0}.
    Ditto for {n1} and the 1-side indices.
    
    The procedure sorts the list of reduced rungs by {msm_rung_lex_compare},
    and discards duplicated rungs.
    
    This procedure made sense when there were circular sequences. 
    Currently it is the same as {msm_pairing_get_rungs}. */

int32_t msm_pairing_count_common_rungs(msm_pairing_t *pa, msm_pairing_t *pb, int32_t n0, int32_t n1);
  /* Given two pairings {pa,pb}, returns the number of rungs that are
    shared between them. 
    
    The rungs considered by this procedure are those returned by
    {msm_pairing_collect_rungs(pa,n0,n1)} and
    {msm_pairing_collect_rungs(pb,n0,n1)}. The result is the number of
    rungs that belong to both sets. */

/* STEPS
  
  A /step/ in a pairing is a pair of consecutive rungs.
  
  A pairing {pr} with {ng} rungs has {ng-1} steps 
  {pr[i] -> pr[i+1]} for {i} in {0..ng-2}. */
  
int32_t msm_pairing_num_steps(msm_pairing_t *pr);
  /* Returns the number steps in the pairing {pr}, or -1 if {pr} is empty. */

/* MAPPING PAIRINGS BETWEEN SCALES */

msm_pairing_t *msm_pairing_map_gen(msm_pairing_t *pr, msm_rung_map_proc_t *map);
  /* Given a pairing {pr}, returns a new pairing {q}
    such that {q[i] = map(pr[i])} for all valid {i}. The result 
    has the same rung count as {pr}. */
    
msm_pairing_t *msm_pairing_map
  ( msm_pairing_t *pr, 
    msm_seq_desc_t *s0_old,
    msm_seq_desc_t *s1_old,
    msm_seq_desc_t *s0_new,
    msm_seq_desc_t *s1_new
  );
  /* Returns a new pairing where each rung {g} of {pr} is mapped by 
    {msm_rung_map(g,s0_old,s1_old,s0_new,s1_new)}. */

/* INCREASING AND ATOMIC PAIRINGS
  
  A pairing is /strictly increasing/ iff every one of
  its steps is strictly increasing.
  
  A strictly increasing pairing can be visualized as a finite
  path in the integer plane, whose steps advance by at least 1 along
  each axis.  */

bool_t msm_pairing_is_increasing(msm_pairing_t *pr, int32_t minIncrEach, int32_t minIncrSum, bool_t atomic, bool_t die);
  /* Returns TRUE iff every step or {pr} advances at least {minIncrEach}
    positions on each side, and a total of at least {minIncrSum} on both
    sequences together. If {atomic} is true, also requires that the
    increment be exactly 1 on at least one side. If these constraints
    are not satisfied, then, if {die} is true, fails with an error
    message, else returns FALSE silently. */

msm_pairing_t *msm_pairing_make_increasing(msm_pairing_t *pr, int32_t minIncrEach, int32_t minIncrSum);
  /* Creates a pairing that contains a subset of the rungs of {pr},
    chosen so that the 0-side and 1-side increments of every step (1) are
    at least {minIncrEach} and (2) add at least {minIncrSum} or more.
    See {msm_rung_vec_make_increasing} for details. */

msm_pairing_t *msm_pairing_interpolate(msm_pairing_t *pr);
  /* Given a strictly increasing pairing {pr}, returns a minimal
    atomic, strictly increasing pairing {q} that contains
    the rungs of {pr} as a subsequence. This is achieved by inserting
    additional rungs between every two successive rungs of {pr}, as
    needed, in the straightest possible way. */

/* PAIRING BETWEEN SEQUENCES 
  
  In this section, each rung of a pairing is understood to be a 
  pair of indices {(i0,i1)} between two abstract sequences {s0,s1},
  respectively with {n0} and {n1} elements, numbered from 0. 
  
  The index {i0} (into the {s0} sequence) must be
  in the range {0..n0-1}, and {i1} must be in {0..n1-1}. */

typedef void msm_pairing_perfect_use_proc_t(int32_t i0, int32_t i1, int64_t ng);
  /* A client procedure that processes a perfect pairing.
    The pairing consists of pairs {(i0+k,i1+k)} for {k} in {0..ng-1}. */

void msm_pairing_enum_alignments
  ( int32_t n0,       /* Length of 0-side sequence. */
    int32_t n1,       /* Length of 1-side sequence. */ 
    msm_pairing_perfect_use_proc_t *use
  );
  /* Enumerates all possible maximal perfect pairings between two
    abstract sequences. The procedure calls {use(i0,i1,ng)} on
    each alignment.
    
    Each pairing will match a prefix of {s0} with a suffix of {s1}, or
    vice-versa, or the whole of {s0} with an internal segment of {s1},
    or vice-versa. */
    
int64_t msm_pairing_perfect_max_length(int32_t n0, int32_t n1);
  /* Maximum length of a perfect pairing between the indices of two 
    abstract sequences; namely, {min(n0,n1)} rungs. The numbers {n0} 
    and {n1} must be positive. */

/* RANDOM PAIRINGS */

msm_pairing_t *msm_pairing_throw
  ( int32_t n0,
    int32_t n1, 
    int32_t len,
    bool_t atomic,
    bool_t diagonal,
    double skipProb
  );
  /* Generates a random pairing {pr}. The 0-side of every rung will be
    restricted to {0..n0-1}, and the 1-side to {0..n1-1}.
    The parameters {n0} and {n1} must be positive. 
    
    The pairing will cover approximately {len} positions on each
    sequence (it may vary somewhat, depending on other constraints). The 
    parameter {len} must not exceed the counts {n0} and {n1}. In any case,
    the number of rungs is chosen by the procedure, normally close to
    {len}.
    
    The pairing will be (strictly) increasing along both sequences.
    If {atomic} is TRUE, the pairing will also be atomic.
    
    If {diagonal} is TRUE, the endpoints are chosen near the diagonal
    of the rectangle {{0..n0-1}×{0..n1-1}}. The parameter {skipProb}
    controls the deviation of the pairing (seen as a path in the XY
    grid) from the straight line connecting the selected end rungs. */

msm_pairing_t *msm_pairing_sub_throw(msm_pairing_t *pr, int32_t len);
  /* Extracts from pairing {*pr} a random set of consecutive rungs that
    span {2*len} positions, total, in the two
    sequences. The total span may be less than {2*len} if {pr} is not atomic. */

/* PAIRING I/O */

void msm_pairing_write(FILE *wr, msm_pairing_t *pr, int32_t ixSize, bool_t writePairing);
  /* Writes the pairing {*pr} to the file {wr}.
    
    The format is "{ng} {i0} {i1} {S}", where {ng} is the number 
    of rungs in {pr} (which must be positive), {(i0,i1)}
    is the first rung of {pr}, and {S} is a compact
    representation of the steps in {pr} (see below).
    The integers are padded to width {ixSize} and separated by whitespace.
    If {writePairing} is FALSE the {S} field is ommited.
    
    If {pr} is an empty or perfect pairing, {S} is a single asterisk "*".
    Otherwise, {S} consists of an initial '|' character, then a sequence
    of {ng-1} steps, as encoded by {msm_rung_step_write}. */

void msm_pairing_write_full(FILE *wr, msm_pairing_t *pr);
  /* Writes the integer pairs of pairing {*pr} into the file {wr}. */
  
msm_pairing_t *msm_pairing_read(FILE *rd);
  /* Reads from the file {rd} a pairing {pr} betqeen two sequences
    {s0} and {s1}, in the format produced by {msm_pairing_write}. */

/* VALIDATION */

bool_t msm_pairing_is_valid(msm_pairing_t *pr, bool_t die);
  /* Runs consistency checks on the pairing {pr}. In particular, requires
    every step {g-->h} of {pr} to be non-decreasing on both sides.
    If {die} is true, and the pairing is not valid, fails with an error
    message. If {die} is FALSE, returns TRUE or FALSE normally without
    any messages. */
  
#endif
