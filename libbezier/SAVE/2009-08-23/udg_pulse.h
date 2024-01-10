/* Univariate polynomial splines defined on irregular dyadic grids */
/* Last edited on 2009-08-22 19:06:39 by stolfi */

#ifndef udg_pulse_H
#define udg_pulse_H

#include <bz_basic.h>

#include <udg_grid.h>
#include <udg_spline.h>

/*  
  UNIDIMENSIONAL CIRCULAR UNIFORM SPLINES
  
  A /unidimensional circular spline/ is a polynomial spline defined on
  some finite unidimensional dyadic grid {G}.
  
  We will denote by {\S_c^g[G]} the space of all unidimensional dyadic
  splines over the grid {G}, whose pieces have maximum degree {g}
  ({\geq 0}) and which are continuous to order {c} ({\geq -1})
  everywhere. 
  
  The {k}th derivative of a spline of order {c} is defined everywhere
  if {k <= c}; for {k > c}, it is defined only inside the cells, not
  at the grid vertices. In particular, a spline with continuity {c =
  -1} is undefined at the grid vertices.
  
  PULSES
  
  We define a /pulse/ as a polynomial spline defined on a regular
  unidimensional grid {G}, infinite or with toroidal (i.e. circular)
  topology.
  
  STANDARD DYADIC PULSE FAMILIES 
  
  We consider here certain /standard dyadic bases/ for the spaces
  {\S_c^g[G]}. Each standard basis is a finite subset of an infinite
  /standard dyadic pulse family/, which is identified by the
  parameters {c} and {g} and by a /kind/ {pkind}.
  
  Every member of such a family is a pulse defined on some level {r}
  of the infinite one-dimensional grid {G*}, namely on a regular grid
  circular with size {gsz = 2^r} for some integer {r\geq 0}. */
  
typedef enum 
  { udg_PK_B = 0, /* Minimal-degree B-spline-like pulses. */
    udg_PK_H = 1, /* Hermite-style pulses. */
    udg_PK_N = 2  /* Compact unit-partition pulses. */
  } udg_pulse_kind_t;

#define udg_pulse_kind_MIN (udg_PK_B)
#define udg_pulse_kind_MAX (udg_PK_N)
  
/* 
  MOTHER PULSES
  
  All pulses of the family are derived from a small number of /mother
  pulses/ through domain scaling, shifting, and reduction to a
  circular grid topology, as explained below. */

typedef uint32_t udg_pulse_mother_count_t; 
  /* A count of mother pulses. */

typedef int32_t udg_pulse_mother_index_t;  
  /* Identifies a mother pulse. */

udg_pulse_mother_count_t udg_pulse_mother_count(udg_pulse_kind_t pkind, udg_cont_t c, udg_degree_t g);
  /* Number of distinct mother pulses in the family defined by pulse
    kind {pkind}, continuity {c}, and degree {g}. Returns zero iff the
    combination {pkind,c,g} is invalid. */

typedef struct udg_pulse_family_t
  { udg_pulse_kind_t pkind;       /* Pulse kind. */
    udg_cont_t c;                 /* Continuity class. */
    udg_degree_t g;               /* Degree. */
    udg_pulse_mother_count_t nmp; /* (READONLY) {= udg_pulse_mother_count(pkind,c,g)}. */
  } udg_pulse_family_t;
  /* Identifies a standard dyadic pulse family from which one can pick a 
    basis for {\S_c^g[G]}, for any finite dyadic grid {G}. */
    
udg_pulse_family_t udg_pulse_family(udg_pulse_kind_t pkind, udg_cont_t c, udg_degree_t g);
  /* Creates a {udg_pulse_family_t} from {c,g,pkind}. The number
    of mother pulses {nmp} is computed from the given parameters. */
     
/*
  RANK, SHIFT, AND INDEX
  
  Assuming that the the root interval {R = [a_b]} is fixed, and a standard
  dyadic pulse family {fam} has been specified, individual elements of that
  family are identified by three additional parameters:

    * a /grid size/ {gsz > 0};
    * a /shift/ {pos} in the range {0..gsz-1};
    * a /mother pulse index/ {pix} in the range {0..fam.nmp-1}.
    
  For dyadic pulses, the grid size is always a power of two, {gsz =
  2^r}; but this interface does not assume it. The grid consists of
  the interval {R} divided into {gsz} equal intervals, each having
  width {delta = (b-a)/gsz}, taken with circular topology.
  
  The shift parameter {pos} is an integer number, to be interpreted modulo
  {gsz}. The /anchor cell/ of the pulse is the cell with position {pos}. */

typedef struct udg_pulse_t 
  { udg_pulse_mother_index_t pix;    /* The mother pulse's index. */
    udg_grid_size_t gsz;             /* Grid size {= 2^r}. */
    udg_grid_pos_t pos;              /* Shift (position of first support interval). */
    udg_grid_size_t msz;             /* (READONLY) Size of mother pulse's support (see below). */
  } udg_pulse_t; 
  /* A {udg_pulse_t} identifies one specific pulse within a pulse family
    (which must be specified separately) on a circular grid with {gsz} cells.
    The pulse is the mother pulse with index {pix}, shifted so that its 
    support begins with cell {pos}, folded-over if needed.
    The field {msz} is a fixed function of {fam} and {pix}. */

udg_pulse_t udg_pulse
  ( udg_pulse_family_t *fam, 
    udg_pulse_mother_index_t pix, 
    udg_grid_size_t gsz, 
    udg_grid_pos_t pos
  );
  /* Returns the pulse from family {fam} with mother pulse index {pix}, grid size
    {gsz}, and shift {pos}. The mother pulse size {msz} is computed
    internally. */

/* 
  STANDARD DYADIC PULSE SUPPORT
  
  The /support/ of a pulse is the shortest list of consecutive grid cells
  that, together with the vertices between them, contains all the
  arguments where the pulse is nonzero. The number of those cells is
  the pulse's /support count/.
  
  For large enough grid size {gsz}, the support count {psz} of a standard
  dyadic pulse is strictly less than the grid size {gsz}, and
  depends only on the family {fam}, and {pix}, not on
  {r} and {pos}. The pulse is nonzero only on the grid cells with
  position {pos..pos+psz-1}, modulo {gsz}.
  
  For small {r}, the count {psz} may be reduced to {gsz} because of
  self-overlap to (see below); in that case the pulse's support is the
  whole grid. 
  
  In any case, the anchor cell of the pulse is the cell that
  corresponds to the lowest cell in the mother pulse's support. */

udg_grid_size_t udg_pulse_supp_count
  ( udg_pulse_family_t *fam,
    udg_pulse_mother_index_t pix, 
    udg_grid_size_t gsz
  );
  /* Returns the number {psz} of consecutive cells in the support of a
    standard dyadic pulse with family {fam}, mother pulse index {pix}, and
    grid size {gsz} such that {gsz = 2^r}. Takes self-overlap into
    account, so that the result never exceeds {gsz}. */

udg_grid_size_t udg_pulse_max_supp_count ( udg_pulse_family_t *fam );
  /* Returns the maximum of {udg_pulse_supp_count(fam,pix,gsz)} over 
    all valid pulse indices {pix} in family {fam} and 
    sufficiently large {gsz}. */

/* 
  STANDARD PULSE EVALUATION
  
  The following procedures evaluate a standard dyadic pulse and all
  its derivatives up to order {ord}, at a specified argument {x}. The
  pulse has family {fam}, mother pulse index {pix}, grid size {gsz} such
  that {gsz=2^r}.

  The {eval} procedures take wrap-around into account when appropriate (see
  below). The result is returned in {f[0..ord]}: the pulse's value in
  {f[0]}, and its {k}th derivative in {f[k]}, for {k=1..ord}. Beware
  that the result may be undefined when {ord > c} and {x} is a
  multiple of {wd/gsz}. */

#define udg_pulse_MAX_DEGREE (7)
  /* A limit on the degree for pulse evaluation, just for safety. */

#define udg_pulse_MAX_DIFF_ORD (udg_pulse_MAX_DEGREE + 1)
  /* A limit on the derivative order, just for safety. */

void udg_pulse_eval
  ( udg_pulse_family_t *fam,
    udg_pulse_t *p,
    interval_t *R,
    double x,
    udg_cont_t ord,
    double *f
  );
  /* Computes the value and derivatives {f[0..ord]} of the pulse {p}
    from family {fam} at the point {x}. Assumes a circular grid with
    size {p.gsz} and root interval {R}, which defaults to {(0_1)}
    if {R == NULL}. */

void udg_pulse_eval_cell_relative
  ( udg_pulse_family_t *fam,
    udg_pulse_mother_index_t pix, 
    udg_grid_size_t gsz,
    double x,
    udg_cont_t ord,
    double *f
  );
  /* Computes the value and derivatives of pulse from family {fam}
    with the given grid size {gsz} and shift {pos = 0}, at argument
    {x}. Assumes a circular grid with step {delta = 1}, i.e. where the
    root interval {R} is {[0 _ gsz)}. */

void udg_pulse_to_bezier
  ( udg_pulse_family_t *fam,
    udg_pulse_t *p,
    udg_grid_pos_t ix,
    double *bz
  );
  /* Extracts the Bézier coefficients of one piece of the standard
    dyadic pulse {p} from family {fam}. Assumes a circular grid with size {gsz}.
    (The result does not depend on the root interval {R}.) The piece is
    the restriction of the spline to the grid cell with position {ix},
    i.e. to the interval {(ix*wd/gsz _ (ix+1)*wd/gsz)}. Takes wrap-around
    into account. The result is nonzero only if {ix - pos} is congruent
    modulo {gsz} to one of the integers {0..psz-1}, where {psz} is the
    pulse's support count. */

/*
  SELF-SIMILARITY AND MOTHER PULSES
  
  Dyadic pulses that differ only in grid size and shift are related to
  each other by argument scaling and translation (like wavelets). The
  other four parameters -- kind {pkind}, continuity {c}, degree {g}, and
  mother pulse index {pix} -- define the basic shape of the pulse, and
  thus will be called the /shape parameters/.
  
  All dyadic pulses with the same shape parameters are derived from a
  single /mother pulse/ {w(x)}. The mother pulse is a polynomial
  spline of continuity order {c} and degree {g}, defined on the
  *non-periodic* integer grid -- the uniform infinite grid on the real
  line with integer vertices.
  
  MOTHER SUPPORT COUNT
  
  The mother pulse {w} has always bounded support, spanning {msz}
  consecutive grid cells (unit intervals) and possibly their inferior
  corner vertices. The number {msz} is the /mother's support count/ for
  pulses of that shape, and is an upper bound for the support size
  {psz} of any derived pulse. By convention, the support of a mother
  spline is always the interval {(0 _ msz)}, that is, cells
  {0..msz-1} of the infinite unit grid.
  
  Mother pulses with continuity order {c = -1} are undefined for
  integer arguments; and for {c \geq 0} the mother pulse's value
  must be zero for arguments {0} and {msz}. */
 
udg_grid_size_t udg_pulse_mother_supp_count
  ( udg_pulse_family_t *fam, 
    udg_pulse_mother_index_t pix
  );
  /* Computes the support count {msz} of the mother spline of 
    family {fam} and mother pulse index {pix}. */

/* MOTHER PULSE EVALUATION

  The following procedures refer to the mother pulse of  
  family {fam}, and mother pulse index {pix}. */

void udg_pulse_mother_eval
  ( udg_pulse_family_t *fam, 
    udg_pulse_mother_index_t pix, 
    double z,
    udg_cont_t ord,
    double *w
  );
  /* Evaluates the mother pulse of family {fam} and mother index {pix},
    and all its derivatives up to order {ord}, at the argument {z}.
    The pulse is specified by its
    
    Returns the result in {w[0..ord]}: the pulse's value in {w[0]},
    and its {k}th derivative in {w[k]}. Beware that {w[k]} may be
    undefined when {k > c} and {z} is an integer. All these values are
    zero if {z} is outside the interval {[0 _ msz]} where 
    {msz = udg_pulse_mother_supp_count(fam,pix)}. */
    
void udg_pulse_mother_to_bezier
  ( udg_pulse_family_t *fam, 
    udg_pulse_mother_index_t pix,
    udg_grid_pos_t ix,
    double *bz
  );
  /* Stores into {bz[0..g]} the Bézier coefficients of the mother
    pulse in the interval {(ix _ ix+1)}. The coefficients will be zero
    if {ix < 0} or {ix >= msz}, where 
    {msz = udg_pulse_mother_supp_count(fam,pix)}. */
 
/* 
  PULSE WRAP-AROUND
  
  Let {w} be a mother pulse with family {fam} and mother index {pix},
  and {p} be a dyadic pulse of family {fam}, defined on a regular
  periodic grid with root interval {R = [a_b)}. The pulse {p} is
  obtained by shifting the mother pulse's domain by {p.pos}, wrapping
  it around the interval {[0 _ p.gsz)} and then changing the
  argument's scale so that {[0 _ p.gsz)} is mapped to {R}. More
  precisely
  
    {p(x) = SUM { w(z) : z \eqv z - p.pos (mod p.gsz) } }  (1)
   
  where {z = p.gsz*(x-a)/(b-a)}. Since there are at most
  {ceil(p.msz/p.gsz)} arguments {z} such that {w(z)} is non-zero, the
  above sum is finite.
  
  In fact, when the grid size is large enough (namely {p.gsz \geq
  msz}), there is at most one {z} in that interval. Then every dyadic
  pulse {p} with those shape parameters looks like a scaled and
  translated version of {w}, and its support consists of {psz =
  min(p.msz,p.gsz)} consecutive cells, whose positions belong to the
  set {p.pos..p.pos+p.psz-1} taken modulo {p.gsz}.
  
  Note also that while the mother pulse {w} is not periodic, the
  derived dyadic pulses are periodic, with fundamental cell {R}.

  When the grid size {p.gsz} is too low (namely {p.gsz < p.msz}), the
  pulse's value may be the sum of {w(z)} for two or more several {z}
  within the mother pulse's support. Intuitively, the mother pulse
  gets folded over itself when the grid is too small; in that case,
  the overlapped parts are added together to give the dyadic pulse
  {p}. */

/* PRINTOUT: */

char udg_pulse_kind_to_char(udg_pulse_kind_t pkind);
  /* Returns a character that identifies kind {pkind}:
    'B' for {udg_PK_B}, 'H' for {udg_PK_H}, etc.)*/

void udg_pulse_family_print(FILE *wr, udg_pulse_family_t *fam);
  /* Prints the pulse family {fam}, in the format
    "{pkind}_{c}^{g}" */

void udg_pulse_print(FILE *wr, udg_pulse_family_t *fam, udg_pulse_t *p);
  /* Prints the pulse {p} from family {fam}, in the format
    "{pkind}_{c}^{g}[{pix}]:{pos}:{msz}/{gsz}", where {gsz} is the
    grid layer size, {pos/gsz} is the lower end of the support,
    {msz/gsz} is the width of the support (before wrap-around).
    The mother index part "[{pix}]" is supressed if there is only one
    mother pulse for that {pkind,c,g} combination. */

#endif
