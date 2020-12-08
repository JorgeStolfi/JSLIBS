#ifndef msm_rung_H
#define msm_rung_H

/* A rung is a pair of integer indices. */
/* Last edited on 2013-10-22 06:52:46 by stolfilocal */

#define msm_rung_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <msm_seq_desc.h>

#include <sign.h>
#include <vec.h>

#include <stdint.h>

/* RUNGS
  
  A /rung/ is a pair of integers. It can be seen as a point of the
  integer plane, so its two components are called the /X coordinate/
  and /Y coordinate/, or /0-coordinate/ and /1-coordinate/.
  
  Each coordinate of a rung is usually meant to be an index into some
  sequence; but this interpretation is not used to this interface,
  exacept as a motivation. */

typedef struct msm_rung_t { int32_t c[2]; } msm_rung_t; 

#define msm_rung_X(g) ((g).c[0])
#define msm_rung_Y(g) ((g).c[1])
  /* The X- and Y-coordinates of a rung. */
  
#define msm_rung_R(g) ((g).c[0]+(g).c[1])
#define msm_rung_S(g) ((g).c[0]-(g).c[1])
  /* The R- and S-coordinates of a rung (along the main
    diagonal and across it, respectively). */

#define msm_rung_none (msm_rung_t){{ INT_MAX, INT_MAX }}
  /* A special {msm_rung_t} value mening /no rung/. */

bool_t msm_rung_is_none(msm_rung_t *g);
  /* TRUE iff {(*g) == msm_rung_none}. */

bool_t msm_rung_equal(msm_rung_t ga, msm_rung_t gb);
  /* TRUE iff {ga == gb}. */

sign_t msm_rung_lex_compare(msm_rung_t ga, msm_rung_t gb);
  /* Compares two rungs {ga,gb} by lexicograpic order. Namely, compares
    the X coordinates, and breaks ties comparing the Y coordinates.
    The result is {NEG} (-1) if {ga} precedes {gb}, {POS} (+1) 
    if {ga} follows {gb}, and {ZER} (0) if {ga == gb}. */
  
sign_t msm_rung_strict_compare(msm_rung_t ga, msm_rung_t gb);
  /* Compares two rungs {ga,gb} by strict domination order. The result
    is {NEG} (-1) if {msm_rung_X(ga) <= msm_rung_X(gb)} and
    {msm_rung_Y(ga) <= msm_rung_Y(gb)} and at least one of the
    inequalities is strict. The result is {POS} (+1) if the same
    condition holds with {ga} and {gb} interchanged. The result is
    {ZER} (0) in all other cases, namely {ga == gb} or the X- and
    Y-orders are inconsistent. */
  
msm_rung_t msm_rung_throw(msm_rung_t glo, msm_rung_t ghi, bool_t diagonal);
  /* Selects a rung at random whose X- and Y-coordinates lie between
    those of {glo} and {ghi}. If {diagonal}, the point will be on the
    approximate diagonal of that rectangular region.
    (TO FIX!!! {diagonal} should be a probability of insert/delete.) */

msm_rung_t msm_rung_throw_middle(msm_rung_t glo, msm_rung_t ghi, double skipProb);
  /* Selects a random rung {g} whose X- and Y-coordinates lie strictly between
    those of {glo} and {ghi}.
    
    The R-coordinate of {g} will be
    {(msm_rung_R(glo)+msm_rung_R(ghi))/2}, randomly rounded. The steps
    {glo-->g} and {g-->ghi} will be (strictly) increasing in both
    coordinates. Therefore, the step {glo->ghi} must increase by at
    least 2 in both coordinates.
    
    The parameter {skipProb} specifies the probability of asymmetric
    steps on a virtual path from {glo} to {ghi}; the rung {g} is taken
    to be the middle rung of that path. (The actual probability may be
    greater than that, depending on {glo,ghi}.) */

vec_typedef(msm_rung_vec_t,msm_rung_vec,msm_rung_t);

msm_rung_vec_t msm_rung_vec_throw(msm_rung_t gini, msm_rung_t gfin, bool_t atomic, double skipProb);
  /* Generates a random sequence of rungs {gv} 
    that starts with rung {gini} and ends with rung {gfin}.
    
    If {gini == gfin}, the vector will have only that rung. Otherwise,
    the step {gini-->gfin} must be (strictly) increasing; and all
    steps of the output rung sequence will be strictly increasing.
    If {atomic} is TRUE, they will be atomic, too.
    
    The parameter {skipProb} controls the deviation of the pairing
    from a straight line. (The actual deviation may be greater than
    that, depending on {gini} and {gfin}.) */

/* STEPS

  A /step/ is a pair of rungs (intended to be between the same two sequences).
  
  A step {g0-->g1} is /atomic/ iff at least one coordinate changes by {±1}.
  
  A step {g0-->g1} is /strictly increasing/ iff both coordinates increase
  by at least 1; that is, {g1.X - g0.X >= 1} and {g1.Y - g0.Y >= 1}.
  
  A step {g0-->g1} is /diagonal/ (a /hop/) if the X and Y coordinates
  change by the same amount; that is, {g1.X - g0.X == g1.Y - g0.Y}.
  A /perfect/ step is a diagonal step with increment {+1}. 
  
  The step is /X-stationary/ if {g0.X == g1.X}, and /Y-stationary/
  if {g.Y == h.Y}.  */
  
bool_t msm_rung_step_is_increasing(msm_rung_t g0, msm_rung_t g1, int minIncrEach, int minIncrSum, bool_t die);
  /* Returns TRUE iff the step {g0-->g1} advances least {minIncrEach}
    positions on each side, and a total of at least {minIncrSum} on both
    sequences together. Othwerwise, if {die} is TRUE, fails with an
    error message; if {die} is FALSE, returns FALSE silently. */

bool_t msm_rung_step_is_atomic(msm_rung_t g0, msm_rung_t g1, bool_t die);
  /* Returns TRUE iff the step {g0-->g1} is atomic.  Othwerwise, if {die} is TRUE,
    fails with an error message; if {die} is FALSE, returns FALSE silently. */

bool_t msm_rung_step_is_perfect(msm_rung_t g0, msm_rung_t g1, bool_t die);
  /* Returns TRUE iff the step {g0-->g1} is perfect.  Othwerwise, if {die} is TRUE,
    fails with an error message; if {die} is FALSE, returns FALSE silently. */

int msm_rung_step_span_increment(msm_rung_t g0, msm_rung_t g1, int j);
  /* If {j} is 0 or 1, returns the signed number of positions spanned
    on side {j} of the rungs by the step from {g0} to {g1}, not
    counting one of the rungs; that is, {g1.c[j] - g0.c[j]}. If {j} is
    {-1}, returns the sum of the spans on both sides; that is,
    {(g1.c[0] - g0.c[0]) + (g1.c[1] - g0.c[1])}. */

sign_t msm_rung_step_break_tie(msm_rung_t a0, msm_rung_t a1, msm_rung_t b0, msm_rung_t b1);
  /* Returns {+1} iff the step {a == a0-->a1} is better than the step
    {b == b0-->b1} according to geometric criteria alone; returns
    {-1} if {b} is better than {a}; and returns 0 if both are
    geometrically equivalent.
    
    This procedure can be used to break the tie when choosing between
    two steps which happen to be equivalent by other (more important)
    criteria.
    
    More precisely, the procedure returns TRUE if the step {a0-->a1}
    is shorter than {b0-->b1} in the L1 metric, or if both have the
    same length but {a0-->a1} is closest to the diagonal. Any steps to or
    from {msm_rung_none} are worse than any steps between two valid
    rungs. */

void msm_rung_interpolate(msm_rung_t g0, msm_rung_t g1, int *ngP, msm_rung_vec_t *gv);
  /* Appends to {gv[0..*ngP-1]} a series of rungs that interpolate
    between rung {g0} (exclusive) and rung {g1} (inclusive), so
    that the steps are atomic, all in the same direction as
    {g0-->g1}, and the resulting path in the X-Y grid is as straight as possible. 
    
    If {g0 == g1}, the procedure is a no-op; otherwise the rungs must
    differ in both coordinates. The procedure increments {*ngP} and
    expands {*gv} as needed. */

msm_rung_vec_t msm_rung_vec_interpolate(msm_rung_vec_t *gv);
  /* Returns a copy of {gv}, with additonal rungs inserted in its steps
    (with {}) so that the result is increasing and atomic.
    
    The steps of {gv} must be strictly increasing or strictly
    decreasing in each coordinate. The new rung vector will begin and
    ednd with the same rungs as {gv}, and will contain all its
    rungs. */

msm_rung_vec_t msm_rung_vec_make_increasing(msm_rung_vec_t *gv, int minIncrEach, int minIncrSum);
  /* Returns a list of rungs that contains a subset of the rungs of {gv},
    chosen so that the X and Y increments of every step (1) are
    at least {minIncrEach} and (2) add at least {minIncrSum} or more.  
    
    The rung vector {gv} must not have any steps with negative increment
    in X or Y. If the total span of {gv} is less than {minIncrSum}, or
    the span on either side is less than {minIncrEach}, the result is a
    single rung, taken somehere in the middle of {gv}. Otherwise, the
    result has the same initial and final rungs as {gv}. */


msm_rung_vec_t msm_rung_vec_join(msm_rung_vec_t *gva, int ja, msm_rung_vec_t *gvb, int jb);
  /* Composes two rung vectors {gva} and {gvb} on their sides {ja} and {jb},
    respectively, producing a single rung vector {gvm}.
    
    Specifically, the procedure finds every pair of rungs {ga} in
    {gva} and {gb} in {gvb} such that coordinate {ja} of {ga} matches
    coordinate {jb} of {gb}. For every such instance, the procedure
    appends to {gvm} a rung {(ia,ib)} where {ia} is coordinate {1-ja}
    of {ga} and {ib} is coordinate {1-jb} of {gb}.
    
    In other words, if each rung vector is interpreted as a relation
    from integers to integers, then the procedure computes the
    composition of relation {gva^{1-2*ja}} with relation
    {gvb^{2*jb-1}}, applied in that order.

    The procedure assumes that {gva,gvb} are strictly increasing, 
    and produces a strictly increasing result. */

typedef msm_rung_t msm_rung_map_proc_t(msm_rung_t g);
  /* A client procedure that applies an arbitrary transformation
    to a rung {g}. */

msm_rung_vec_t msm_rung_vec_map_gen(msm_rung_vec_t *gv, msm_rung_map_proc_t *map);
  /* Given a rung list {gv}, returns a new rung list {hv}
    such that {hv[i] = map(gv[i])} for all valid {i}. */

msm_rung_vec_t msm_rung_vec_map
  ( msm_rung_vec_t *gv, 
    msm_seq_desc_t *s0_old,
    msm_seq_desc_t *s1_old,
    msm_seq_desc_t *s0_new,
    msm_seq_desc_t *s1_new
  );
  /* Assumes that the rungs in {gv} relate indices in {s0_old} and in
    {s1_old}. Maps every rung so that it relates the corresponding
    indices of {s0_new} and {s1_new}. Namely, maps the X-coordinate {i0}
    of every rung with {msm_seq_desc_map_index(i0,s0_old,s0_new)}, and its
    Y-coordinate by {msm_seq_desc_map_index(i1,s1_old,s1_new)}.
    
    The sequences {s0_old} and {s0_new} must have the same {id} and {name},
    and ditto for {s1_old} and {s1_new}.
    
    When the mapping implies a coarsening (i.e. when {s0_new.estep <
    s0_old.estep}) the mapping would usually result in fractional rung
    indices. Those are rounded to the nearest integer, with ties rounded
    to even. Note that the rung vector may become non-increasing in that
    case. */

/* RUNG/STEP SCORING */

typedef double msm_rung_score_proc_t
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    msm_rung_t *g
  );
  /* A procedure that computes a numerical ``goodness'' score for a
    single rung {g} between sequences {seq0} and {seq1}. The two sequences
    should usually have the same filtering and {estep}. 
    If {g} is {msm_rung_none}, the function should return 0. */

typedef double msm_rung_step_score_proc_t
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    msm_rung_t *g0,
    msm_rung_t *g1
  );
  /* A procedure that computes a numerical ``goodness'' score for a
    single step between sequences {seq0} and {seq1}, from rung {g0} to
    rung {g1}. The two sequences should usually have the same
    filtering level and {estep}. 
    
    If one of the rungs is {msm_rung_none}, the procedure should
    return some partial score for the partial step. If both are
    {msm_rung_none}, the function should return 0. */

/* RUNG/STEP I/O */

void msm_rung_step_write(FILE *wr, int32_t d0, int32_t d1, int32_t rep);
  /* Writes a compact representation of {rep} steps (and {rep} rungs),
    consisting of an arbitrary step from some unspecified previous rung
    followed by {rep-1} perfect steps. The {rep} argument must be
    positive, and the arguments {d0} and {d1} are the increments in each
    coordinate of the first step (which must not be negative).
    
    The output has the format "({d0},{d1})|{rep}", or just 
    "({d0},{d1})|" if {rep} is 1. The following abbreviations
    are also used to indicate certain initial steps, where
    {D} is a non-negative integer (number of samples skipped
    on one or both sides):
    
      "|"       means "(1,1)|" (a perfect step, skips nothing)
      "/"       means "(1,0)|" (advances 1 in X, repeats Y, skips nothing)
      "\"       means "(0,1)|" (repeats X, advances 1 in Y, skips nothing)
      
      ":|"      means "(2,2)|" (skips 1 sample on both sides)
      "::|"     means "(3,3)|", etc. 
      ":{D}|"   means "({D+1},{D+1})|" (skips {D} samples on both sides)
      
      "'|"      means "(2,1)|" (skips 1 sample in X, advances by 1 in Y)
      "''|"     means "(3,1)|", etc.
      "'{D}|"   means "({D+1},1)|" (skips {D} samples on X, advances 1 in Y)
      
      "'/"      means "(2,0)|" (skips 1 sample in X, repeats Y)
      "''/"     means "(3,0)|", etc.
      "'{D}/"   means "({D+1},0)|" (skips {D} samples on X, repeats Y)
     
      ".|"      means "(1,2)|"  (advances by 1 in X, skips 1 sample in Y)
      "..|"     means "(1,3)|", etc     
      ".{D}|"   means "(1,{D+1})|"(advances by 1 in X, skips {D} samples on Y)
       
      ".\"      means "(0,2)|" (repeats X, skips 1 sample in Y)
      "..\"     means "(0,3)|", etc.
      ".{D}\"   means "(0,{D+1})|" (repeats X, skips {D} samples on Y)
      
    The format changed on 2013-10-16 in an incompatible way.
    Formerly, there was no {rep} count after the last character, the 
    pairs "({d0},{d1}) were written without parentheses,
    and the "'{D}/" or ".{D}\" cases did not have the '\'' 
    or '.' characters. */

void msm_rung_step_read(FILE *rd, int32_t *d0P, int32_t *d1P, int32_t *repP);
  /* Reads from the file {rd} an arbitrary step followed by
    zero or more perfect steps, as encoded by {msm_pairing_step_write}.
    The increments are returned in {*d0P,*d1P}, and the total number
    of steps (and rungs) is returned in {*repP}.
    
    The increments {*d0P,*d1P} will be non-negative and the {*repP} count
    will be positive.  Fails noisily if {rd} is at end-of-file or the
    format is invalid. */

#endif

