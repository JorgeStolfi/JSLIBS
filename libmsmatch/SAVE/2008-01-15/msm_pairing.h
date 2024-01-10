#ifndef msm_pairing_H
#define msm_pairing_H

/* Pairings between elements of two open or circular index ranges */
/* Last edited on 2008-01-12 00:36:57 by stolfi */

#define msm_pairing_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#include <msm_basic.h>
#include <msm_rung.h>

#include <sign.h>
#include <vec.h>

#include <stdint.h>

/* PAIRINGS

  A /pairing/ is a sequence of one or more rungs {p[i]}, for some set
  of integer /rung indices/ {i}.
  
  There are two basic kinds of pairings: /open/ and /circular/. 
  
  If {p} is an open pairing, then {p[i]} is defined only for rung
  indices {i} in {0..ng-1}, for some integer {ng >= 0}.
  
  If {p} is a circular pairing, then {p[i]} is defined for any integer
  rung index {i}. However, the rungs are /integral-periodic/, meaning
  that there is a positive integer {ng} (the /rung index period/) and
  integers {DX,DY} (the /coordinate periods/) such that {p[i+ng] =
  p[i]+(DX,DY)} for every {i}. The /fundamental rungs/ of a circular
  pairing are the rungs with indices {0..ng-1]}. */
  
typedef struct msm_pairing_t msm_pairing_t;

bool_t msm_pairing_is_circular(msm_pairing_t *p);
  /* Returns TRUE iff {p} is a circular pairing. */

int msm_pairing_num_rungs(msm_pairing_t *p);
  /* Returns the number {ng} of rungs in the pairing {p}. 
    If {p} is circular, the result is the rung index period. */

msm_rung_t msm_pairing_period(msm_pairing_t *p);
  /* Returns the coordinate periods {PX,PY} of {p}, packed as a rung.
    If {p} is open, the result is {(0,0)}. */

msm_rung_t msm_pairing_get_rung(msm_pairing_t *p, int i);
  /* Obtains the rung with index {i} in the pairing {*p}. If {p} is
    open, the index {i} must be in the range {0..ng-1}, where {ng} is
    the number of rungs. If {p} is circular, the index {i} can be any
    integer. */

bool_t msm_pairing_equal(msm_pairing_t *pa, msm_pairing_t *pb, bool_t die);
  /* Returns TRUE if the pairings {pa} and {pb} are equal.
    Otherwise, if {die} is true, fails with an error message; if {die}
    is FALSE, returns FALSE silently. */

msm_rung_vec_t msm_pairing_get_rungs(msm_pairing_t *p);
  /* Returns a new rung vector with a copy of all the rungs of {p}.
    
    If {p} is open, the result will contain only the {ng} fundamental
    rungs of {p}, where {ng = msm_pairing_num_rungs(p)}.
    
    If {p} is circular, the result will have {ng+1} rungs; elements
    {0..ng-1} are the fundamental rungs of {p}, starting with an
    arbitrary rung {gIni} in the and element {ng} is {gIni+(DX,DY)},
    where {(DX,DY)} is the coordinate period vector of {p}. */

msm_pairing_t *msm_pairing_from_rung_vec(msm_rung_vec_t *gv, bool_t circp);
  /* Creates a pairing from a rung vector {gv}.
    
    If {circp} is false, the pairing will be open, and it will have
    {nv = gv.nel} rungs, namely {gv.el[0..nv-1]}.
    
    If {circp} is true, the pairing will be circular. Its rung index
    period will be {ng = gv.nel-1}, the rungs with indices {0..ng-1}
    will be {gv.el[0..ng-1]}, and the coordinate period vector
    will be {gv[ng] - gv[0]}.
    
    In either case, the order and multiplicity of the rungs is
    preserved. The pairing will share the vector {gv->el}. */

void msm_pairing_free(msm_pairing_t *p);
  /* Frees all data records associated with {p},
    including {(*p)} and its rung vector, if any. */

/* PERFECT PAIRINGS
  
   A pairing is /perfect/ if the difference between any two consecutive 
   rungs is {(+1,+1)}. */

bool_t msm_pairing_is_perfect(msm_pairing_t *p, bool_t die);
  /* TRUE iff {p} is a perfect pairing. Otherwise, if {die} is true,
    fails with an error message, else returns FALSE silently. */

msm_pairing_t *msm_pairing_perfect(msm_rung_t *gini, int ng, bool_t circp);
  /* Creates a perfect pairing with rungs {gini + (k,k)} for {k} in {0..ng-1}.
    
    If {circp} is false, the pairing will be an open one, with those
    {ng} rungs only. If {circp} is true, the pairing will be circular,
    with rung index period {ng} and coordinate periods {(ng,ng)}. */

msm_rung_vec_t msm_pairing_collect_rungs(msm_pairing_t *p, int nx, int ny);
  /* Returns the set of all distinct fundamental rungs of {p},
    namely rungs {p[0..ng-1]} where {ng = msm_pairing_num_rungs(p)},
    reduced modulo {nx,ny}.
    
    The parameters {nx} and {ny} must be non-negative. If {nx} is
    positive, the X coordinates of every rung are reduced modulo {nx}.
    (If {p} is circular, this reduction makes sense only if its X period {DX}
    is a multiple of {nx}.)  Ditto for {ny} and the Y coordinates.
    
    The procedure sorts the list of reduced rungs by {msm_rung_lex_compare},
    and discards duplicated rungs. */

int msm_pairing_count_common_rungs(msm_pairing_t *pa, msm_pairing_t *pb, int nx, int ny);
  /* Given two pairings {pa,pb}, returns the number of rungs that are
    shared between them. 
    
    The rungs considered by this procedure are those returned by
    {msm_pairing_collect_rungs(pa,nx,ny)} and
    {msm_pairing_collect_rungs(pb,nx,ny)}. The result is the number of
    rungs that belong to both sets. */

/* STEPS
  
  A /step/ in a pairing is a pair of consecutive rungs.
  
  An open pairing {p} with {ng} rungs has {ng-1} steps 
  {p[i] -> p[i+1]} for {i} in {0..ng-2}.
  
  A circular pairing {p} with rung index period {ng} has infintely
  many steps, but exactly {ng} /fundamental steps/ {p[i] -> p[i+1]}
  for {i} in {0..ng-1}. */
  
int msm_pairing_num_steps(msm_pairing_t *p);
  /* Returns the number steps in the pairing {p}. If {p} is circular,
    the result is the number of fundamental steps. */

/* MAPPING PAIRINGS BETWEEN SCALES */

msm_pairing_t *msm_pairing_map_coordinates(msm_pairing_t *p, msm_rung_map_proc_t *map);
  /* Given a pairing {p}, returns a new pairing {q}
    such that {q[i] = map(p[i])} for all valid {i}.
    
    If {p} is open, {q} is open and has the same rung count. If {p} is
    circular, {q} is circular and has the same rung index period. */
    
msm_pairing_t *msm_pairing_map_to_finer
  ( msm_pairing_t *p, 
    int nxold, 
    int nyold, 
    int nxnew, 
    int nynew, 
    int nwtb
  );
  /* Returns a new pairing whose rung list is 
    {msm_filter_map_rung_vec_to_finer(p->gv,nxold,nyold,nxnew,nynew)}. */

msm_pairing_t *msm_pairing_map_to_coarser
  ( msm_pairing_t *p, 
    int nxold, 
    int nyold, 
    int nxnew, 
    int nynew, 
    int nwtb
  );
  /* Returns a new pairing whose rung list is 
    {msm_filter_map_rung_vec_to_coarser(p->gv,nxold,nyold,nxnew,nynew)}. */

/* MONOTONIC AND UNITARY PAIRINGS
  
  A pairing, open or circular, is /{M}-increasing/ iff every one of
  its steps is {M}-increasing.
  
  An {M}-increasing pairing can be visualized as a finite or infinite
  path in the integer plane, whose steps advance by at least {M} along
  each axis.
  
  Note that a circular pairing is {M}-increasing iff its fundamental steps 
  are {M}-increasing.  */

bool_t msm_pairing_is_increasing(msm_pairing_t *p, bool_t die);
  /* Returns TRUE iff {p} is (strictly) increasing in both sequences. Otherwise,
    if {die} is true, fails with an error message, else
    returns FALSE silently. */

bool_t msm_pairing_is_atomic_increasing(msm_pairing_t *p, bool_t die);
  /* Returns TRUE iff all steps in pairing {p} are (strictly)
    increasing in both sequences, and increase by exactly 1 along one
    or both of them. Otherwise, if {die} is true, fails with an error
    message, else returns FALSE silently. */

msm_pairing_t *msm_pairing_interpolate(msm_pairing_t *p);
  /* Given a 1-increasing pairing {p}, returns a minimal
    1-atomic, 1-increasing pairing {q} that contains
    the rungs of {p} as a subsequence. This is achieved by inserting
    additional rungs between every two successive rungs of {p}, as
    needed, in the straightest possible way.
    
    The result will be circular iff {p} is circular. In that case, the
    result will have the same coordinate period vector as {p}. */

/* PAIRING BETWEEN SEQUENCES 
  
  In this section, each rung of a pairing is understood to be a 
  pair of indices {(ix,iy)} between two abstract sequences {x,y},
  respectively with {nx} and {ny} elements, numbered from 0. 
  
  Each sequence may be either open or circular, as indicated by flags
  {circx} and {circy}. If {circx} is FALSE, the X-coordinate of each
  rung (an index into the {x} sequence) must bein the range {0..nx-1}.
  If {circx} is TRUE, the X-coordinate can be any integer {ix}, which
  is assumed to mean the element with index {imod(ix,nx)} of the {x}
  sequence. The parameters {circy} and {ny} have the same meaning for
  the {Y}-coordinate. */

typedef void msm_pairing_perfect_use_proc_t(int ix, int iy, int64_t ng, bool_t circp);
  /* A client procedure that processes a perfect pairing.
    The pairing consists of pairs {(ix+k,iy+k)} for {k} in {0..ng-1}.
    If {circp} is TRUE, the pairing itself should be assumed to be circular. */

void msm_pairing_enum_alignments
  ( int nx,       /* Number of valid X-coordinates. */
    bool_t circx, /* TRUE if X-range is circular. */
    int ny,       /* Number of valid Y-coordinates. */ 
    bool_t circy, /* TRUE if Y-range is circular. */ 
    msm_pairing_perfect_use_proc_t *use
  );
  /* Enumerates all possible maximal perfect pairings between two
    abstract sequences. The procedure calls
    {use(ix,nx,iy,ny,ng,circp)} on each alignment.
    
    If both sequences are circular ({circx==circy==TRUE}), every
    reported pairing will be circular, with length {lcm(nx,ny)}. If
    only {x} is circular, each pairing will cover {y} completely, and
    run {nx/ny} times around {x}. The analogous statement holds if
    only {y} is circular. If neither sequence is circular, each
    pairing will match a prefix of {x} with a suffix of {y}, or
    vice-versa, or the whole of {x} with an internal segment of {y},
    or vice-versa. */
    
int64_t msm_pairing_perfect_max_length(int nx, bool_t circx, int ny, bool_t circy);
  /* Maximum length of a perfect pairing between the indices of two 
    open or circular abstract sequences.
    
    If {circy} is FALSE or {circy} is FAlSE, all perfect pairings
    between the two coordinate ranges are open. If {circx} and {circy}
    are both FALSE, the longest perfect pairing has {min(nx,ny)}
    fundamental rungs. If only {circx} is FALSE, the result is {nx};
    if only {circy} is FALSE, the result is {ny}.
    
    If both {circx} and {circy} are TRUE, one can have infinite (open)
    perfect pairings. The `longest' perfect pairings are then assumed
    to be the circular perfect pairings whose fundamental rungs are
    all distinct modulo {(nx,ny)}. Those pairings all have
    {lcm(nx,ny)} fundamental rungs.
    
    The numbers {nx} and {ny} must be positive. */

/* RANDOM PAIRINGS */

msm_pairing_t *msm_pairing_throw
  ( int nx,
    bool_t circx, 
    int ny, 
    bool_t circy, 
    int len,
    bool_t circp, 
    bool_t atomic,
    bool_t diagonal,
    double skipProb
  );
  /* Generates a random pairing {p}, starting with a rung in
    the set {{0..nx-1}×{0..ny-1}}. The parameters {nx} and {ny} must
    be positive.
    
    If {circx} is FALSE, the X coordinates of all rungs will be
    restricted to {0..nx-1}; if {circy} is FALSE, the Y coordinates
    are restricted to {0..ny-1}.  
    
    The pairing will cover approximately {len} bases on each sequence
    (it may vary somewhat, depending on other constraints). If either
    sequence is open, {len} must not exceed the corresponding sample
    count ({na} or {nb}). In any case, the number of rungs is chosen
    by the procedure, normally close to {len}.
    
    The pairing will be circular if and only if {circp} is TRUE. In
    that case, at least one of {circx} and {circy} must be TRUE. If
    {circx} is TRUE, the period {PX} will be some multiple of {nx}
    (possibly 0), otherwise it will be 0. Ditto for {circy,PY,ny}.
    
    The pairing will be (strictly) increasing along both sequences.
    If {atomic} is TRUE, the pairing will also be atomic.
    
    If {diagonal} is TRUE, the endpoints are chosen near the diagonal
    of the rectangle {{0..nx-1}×{0..ny-1}}. The parameter {skipProb}
    controls the deviation of the pairing (seen as a path in the XY
    grid) from the straight line connecting the selected end rungs. */

/* PAIRING I/O */

void msm_pairing_write(FILE *wr, msm_pairing_t *p, int ixsize);
  /* Writes the pairing {*p} to the file {wr}.
    
    The format is the number of pairs {ng} (which must be positive),
    the initial rung {ix,iy}, the coordinate period vector {PX,PY} (or
    0,0 if the pairing is open), the final rung {fx,fy}, and the rungs
    {p[0..ng-1]}. The integers are separated by whitespace.
    
    A perfect pairing is encoded by a single '*' character. Other
    pairings are encoded by an initial '|' character, then a sequence
    of {ng-1} steps, as described under {msm_pairing_step_write}. */

msm_pairing_t *msm_pairing_read(FILE *rd);
  /* Reads from the file {rd} a pairing {p} betqeen two sequences
    {x} and {y}, in the format produced by {msm_pairing_write}. */

/* VALIDATION */

bool_t msm_pairing_is_valid (msm_pairing_t *p, bool_t die);
  /* Runs consistency checks on the pairing {p}.
    If {die} is true, and the pairing is not valid,
    fails with an error message.  If {die} is FALSE,
    returns TRUE or FALSE normally without any messages. */

#endif
