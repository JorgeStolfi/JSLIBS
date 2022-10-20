/* Univariate polynomial splines defined on regular circular grids */
/* Last edited on 2011-09-19 01:03:50 by stolfilocal */

#ifndef psp_pulse_H
#define psp_pulse_H

#include <psp_basic.h>
#include <stdint.h>
#include <psp_grid.h>

/*
  PULSES
  
  We define a /pulse/ as a finite element (a polynomial spline with
  finite support) on a regular circular unidimensional grid {G}.
  
  STANDARD PULSE FAMILIES
  
  We consider here certain /standard pulses families/, that can be
  used to build bases bases for spline spaces on regular circular
  unidimensional grids. All the pulses in a standard pulse family are
  obtained from a finite list of /mother pulses/ through domain
  scaling, shifting, and reduction to a circular grid topology, as
  explained below.
  
  A standard pulse family is identified by a
  continuity {c}, a degree {g} and by a /pulse kind/ {pkind}. */
  
typedef enum 
  { psp_PK_B = 0, /* Minimal-degree B-spline-like pulses. */
    psp_PK_H = 1, /* Hermite-style pulses. */
    psp_PK_N = 2  /* Compact unit-partition pulses. */
  } psp_pulse_kind_t;

#define psp_pulse_kind_MIN (psp_PK_B)
#define psp_pulse_kind_MAX (psp_PK_N)
  
typedef uint32_t psp_pulse_mother_count_t; 
  /* A count of mother pulses. */

typedef int32_t psp_pulse_mother_index_t;  
  /* Identifies a mother pulse within a list of mother pulses. */

typedef struct psp_pulse_family_t
  { psp_pulse_kind_t pkind;       /* Pulse kind. */
    psp_cont_t c;                 /* Continuity class. */
    psp_degree_t g;               /* Degree. */
    psp_pulse_mother_count_t nmp; /* (READONLY) {= psp_pulse_mother_count(pkind,c,g)}. */
  } psp_pulse_family_t;
  /* A record that describes a standard pulse family. */
    
psp_pulse_family_t psp_pulse_family(psp_pulse_kind_t pkind, psp_cont_t c, psp_degree_t g);
  /* Creates a {psp_pulse_family_t} from {c,g,pkind}. The number
    of mother pulses {nmp} is computed from the given parameters. */

#define psp_pulse_MAX_DEGREE (7)
  /* A limit on the degree for pulse evaluation, just for safety. */

/*
  INDIVIDUAL STANDARD PULSES
  
  Assuming that the the root interval {R = [a_b]} is fixed, and a standard
  pulse family {fam} has been specified, individual members of that
  family are identified by three additional parameters:

    * a /grid size/ {gsz > 0};
    * a /shift/ {pos} in the range {0..gsz-1};
    * a /mother pulse index/ {pix} in the range {0..fam.nmp-1}.
    
  The grid consists of the interval {R} divided into {gsz} equal
  intervals, each having width {delta = (b-a)/gsz}, taken with
  circular topology.
  
  The shift parameter {pos} is an integer number, to be interpreted modulo
  {gsz}. The /anchor cell/ of the pulse is the cell with position {pos}. */

typedef struct psp_pulse_t 
  { psp_grid_size_t gsz;             /* Grid size. */
    psp_grid_pos_t pos;              /* Shift (position of first support interval). */
    psp_pulse_mother_index_t pix;    /* The mother pulse's index. */
    psp_grid_size_t msz;             /* (READONLY) Size of mother pulse's support (see below). */
  } psp_pulse_t; 
  /* A {psp_pulse_t} identifies one specific pulse within a standard pulse family
    (which must be specified separately) on a circular grid with {gsz} cells.
    The pulse is the mother pulse with index {pix}, shifted so that its 
    support begins with cell {pos}, folded-over if needed.
    The field {msz} is a fixed function of {fam} and {pix}. */

psp_pulse_t psp_pulse
  ( psp_pulse_family_t *fam, 
    psp_pulse_mother_index_t pix, 
    psp_grid_size_t gsz, 
    psp_grid_pos_t pos
  );
  /* Returns the pulse from family {fam} with mother pulse index {pix}, grid size
    {gsz}, and shift {pos}. The mother pulse size {msz} is computed
    internally. */

/* 
  EVALUATION OF STANDARD PULSES
  
  The following procedures evaluate a standard pulse {*p} and all
  its derivatives up to order {ord}, at any specified argument {x}.

  The {eval} procedures take self-overlap into account when appropriate (see
  below). The result is returned in {f[0..ord]}: the pulse's value in
  {f[0]}, and its {k}th derivative in {f[k]}, for {k=1..ord}. Beware
  that the result may be undefined when {ord > c} and {x} is a
  multiple of {wd/gsz}. */

#define psp_pulse_MAX_DIFF_ORD (psp_pulse_MAX_DEGREE + 1)
  /* A limit on the derivative order, just for safety. */

void psp_pulse_eval
  ( psp_pulse_family_t *fam,
    psp_pulse_t *p,
    interval_t *R,
    double x,
    psp_cont_t ord,
    double *f
  );
  /* Computes the value and derivatives {f[0..ord]} of the pulse {p}
    from family {fam} at the point {x}. Assumes a circular grid with
    size {p.gsz} and root interval {R}, which defaults to {(0_1)}
    if {R == NULL}. */

void psp_pulse_eval_cell_relative
  ( psp_pulse_family_t *fam,
    psp_pulse_mother_index_t pix, 
    psp_grid_size_t gsz,
    double x,
    psp_cont_t ord,
    double *f
  );
  /* Computes the value and derivatives of pulse from family {fam}
    with the given grid size {gsz} and shift {pos = 0}, at argument
    {x}. Assumes a circular grid with step {delta = 1}, i.e. where the
    root interval {R} is {[0 _ gsz)}. */

void psp_pulse_to_bezier
  ( psp_pulse_family_t *fam,
    psp_pulse_t *p,
    psp_grid_pos_t ix,
    double *bz
  );
  /* Extracts the Bézier coefficients of one piece of the standard
    pulse {p} from family {fam}. Assumes a circular grid with size
    {gsz}. (The result does not depend on the root interval {R}.) The
    piece is the restriction of the spline to the grid cell with
    position {ix}, i.e. to the interval {(ix*wd/gsz _ (ix+1)*wd/gsz)}.
    Takes wrap-around and self-overlap into account. The result is
    nonzero only if {ix - pos} is congruent modulo {gsz} to one of the
    integers {0..psz-1}, where {psz} is the pulse's support count. */
 
/* 
  PULSE WRAP-AROUND AND SELF-OVERLAP
  
  Let {p} be a pulse with mother pulse {w} and root interval {R =
  [a_b)}. The pulse {p} is obtained by shifting the mother pulse's
  domain by {p.pos}, wrapping it around the interval {[0 _ p.gsz)} and
  then changing the argument's scale so that {[0 _ p.gsz)} is mapped
  to {R}.

  When the grid size {p.gsz} is smaller than the mother pulse's
  support count {p.msz}, the wrap-around causes the mother pulse to
  overlap itself. in that case, by definition, the overlapped parts
  are added together to give the pulse {p}. In other words, the
  pulse's value {p(x)} is the sum of {w(z)} for all {z} that
  correspond to {x} by the domain shifting, fold-over, and scaling.
  More precisely,
  
    {p(x) = SUM { w(z) : z \eqv z0 - p.pos (mod p.gsz) } }  (1)
   
  where {z0 = p.gsz*(x-a)/(b-a)}. Since there are at most
  {ceil(p.msz/p.gsz)} such arguments {z} such that {w(z)} is non-zero,
  the above sum is finite.
  
  In fact, when the grid size is large enough (namely {p.gsz \geq
  p.msz}), there is at most one {z} in that interval. Then the pulse {p}
  looks like the mother pulse {w}, except by somain scaling,
  translation, and wrap-around.
  
  Note also that while the mother pulse {w} is not periodic, the
  derived standard pulses are periodic, with fundamental cell {R}.
  
  STANDARD PULSE SUPPORT
  
  If a pulse {p} has {p.gsz \geq p.msz}, the definitions above imply that its
  support count {psz} is the same as its mother's, namely {p.msz};
  which depends only on the family {fam} and on the mother
  index {p.pix}, not on {p.gsz} and {p.pos}.  In that case, the support
  consists of the cells with positions {p.pos..p.pos+p.msz-1}, modulo {p.gsz}.
  
  If {p.gsz < p.msz}, however, the wrap-around causes the support
  {p.psz} to be less than {p.msz}. By convention, in that case we
  define the pulse's support as being the entire grid, and the suport
  count to be {p.gsz}. (In theory, for an arbitrary mother pulse, the
  wrap-around could lead to complete cancellation in some cells of the
  grid. However, we conjecture that such thing cannot happen if the
  mother pulses generate a minimal support basis.)
  
  In any case, the anchor cell of the pulse is the cell that
  corresponds to the first cell (the interval {[0 _ 1]}) in the mother
  pulse's support. */

psp_grid_size_t psp_pulse_supp_count
  ( psp_pulse_family_t *fam,
    psp_pulse_mother_index_t pix, 
    psp_grid_size_t gsz
  );
  /* Returns the number {psz} of consecutive cells in the support of a
    standard pulse with family {fam}, mother pulse index {pix}, and
    grid size {gsz}. Takes wrap-around into account, so that the
    result never exceeds {gsz}. */

psp_grid_size_t psp_pulse_max_supp_count ( psp_pulse_family_t *fam );
  /* Returns the maximum of {psp_pulse_supp_count(fam,pix,gsz)} over 
    all valid mother pulse indices {pix} in family {fam} and 
    sufficiently large {gsz}. */

/* PRINTOUT: */

char psp_pulse_kind_to_char(psp_pulse_kind_t pkind);
  /* Returns a character that identifies kind {pkind}:
    'B' for {psp_PK_B}, 'H' for {psp_PK_H}, etc.)*/

void psp_pulse_family_print(FILE *wr, psp_pulse_family_t *fam);
  /* Prints the pulse family {fam}, in the format
    "{pkind}_{c}^{g}" */

void psp_pulse_print(FILE *wr, psp_pulse_family_t *fam, psp_pulse_t *p);
  /* Prints the pulse {p} from family {fam}, in the format
    "{pkind}_{c}^{g}[{pix}]:{pos}:{msz}/{gsz}", where {gsz} is the
    grid layer size, {pos/gsz} is the lower end of the support,
    {msz/gsz} is the width of the support (before wrap-around).
    The mother index part "[{pix}]" is supressed if there is only one
    mother pulse for that {pkind,c,g} combination. */

#endif
