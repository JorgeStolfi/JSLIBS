#ifndef box_H
#define box_H

/*  Axis-aligned boxes in {R^d}.  */
/*  Last edited on 2023-02-18 22:14:25 by stolfi */

/* We need to set these in order to get {asinh}. What a crock... */
#undef __STRICT_ANSI__
#define _ISOC99_SOURCE 1
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <set32.h>
#include <interval.h>

/*
  An /{m}-dimensional box of {R^d}/ is the Cartesian product {B} of
  {d} real sets {B_i}, {i = 0,..d-1}; of which {m} are open intervals,
  and {d-m} are singleton sets.
  
  A box {B} of {R^d} is here represented by an array of {interval_ts}
  {B[0..d-1]}, where {B[i]} is interpreted as open if {LO(B[i] < HI(B[i])},
  and as closed (i.e., singleton) if {LO(B[i]) == HI(B[i])}.
  
  A box {B} is empty if any of the intervals {B[0..d-1]} is empty;
  that is, if {LO(B[i]) > HI(B[i])}.
  
  There is only one zero-dimensional box (with {d=0}), and it is NOT empty. */
  
/* DATA TYPES */

typedef uint8_t box_dim_t; 
  /* The dimension of a space, set, basis, etc.. */

#define box_MAX_DIM (30)
  /* Maximum dimension of a box. Must not exceed {set32_ELEM_MAX+1}. */

typedef int8_t box_axis_t; 
  /* Identifies a coordinate axis. We assume that the axes
    of {R^d} are numbered from {0} to {d-1}. */

#define box_NO_AXIS (set32_NO_ELEM)
  /* An invalid {box_axis_t} value. */

/* BOX GEOMETRIC PARAMETERS */

bool_t box_is_empty(box_dim_t d, interval_t B[]);
  /* True iff the box {B0..d-1]} is empty; that is, any of 
    the intervals {B[i]} is empty.  Returns {FALSE} if {d} is zero. */

void box_lo_corner(box_dim_t d, interval_t B[], double p[]);
  /* Stores in {p[0..d-1]} the /inferior corner/ of the box {B[0..d-1]},
    i.e. sets {p[i] = lo(B[i])} for all {i}.  Fails if {B} is empty. */

void box_hi_corner(box_dim_t d, interval_t B[], double p[]);
  /* Stores in {p[0..d-1] the /superior corner/ of box {B[0..d-1]}.
    Fails if {B} is empty. */
  
void box_corner(box_dim_t d, interval_t B[], interval_side_t dir[], double p[]);
  /* Stores in {p[0..d-1]} the corner of the box {B[0..d-1]} which,
    along each axis {i}, lies in the direction {dir[i]} -- namely,
    {p[i] = B[i].end[dir[i]]} for all {i}.
    Fails if {B} is empty. */

void box_center(box_dim_t d, interval_t B[], double p[]);
  /* Stores in {p[0..d-1]} the center of the box {B[0..d-1]} -- namely,
    {p[i] = interval_mid(B[i])} for all {i}.
    Fails if {B} is empty. */

void box_half_widths(box_dim_t d, interval_t B[], double h[]);
  /* Stores in {h[0..d-1]} the half-extent of the box {B[0..d-1]} along each 
    axis {i}.  Namely, if {B} is not empty, it sets {h[i] = interval_rad(B[i])}
    for all {i}.
    
    If {B} is empty, sets all {h[i]} to zero. Note that this is
    inconsistent with {interval_rad}, that would return {-INF}. */

void box_widths(box_dim_t d, interval_t B[], double w[]);
  /* Stores in {w[0..d-1]} the extent of the box {B[0..d-1]} along each 
    axis {i} -- namely, {w[i] = interval(B[i].end[1] - B[i].end[1])/2} for all {i}.
    If {B} is empty, sets all {h[i]} to zeros. */

double box_max_width(box_dim_t d, interval_t B[]);
  /* Returns the maximum extent of box {B[0..d-1]} along any axis.
    If {B} is empty or {d} is zero, returns zero. */

double box_radius(box_dim_t d, interval_t B[]);
  /* Euclidean radius of box {B[0..d-1]}. If {B} is empty or {d} is zero,
    returns zero. */

bool_t box_equal(box_dim_t d, interval_t A[], interval_t B[]);
  /* Returns true iff {A} and {B} are both empty, or both non-empty
    and identical. */

/* BOX CREATION AND MODIFICATION */

void box_empty(box_dim_t d,  interval_t C[]);
  /* Sets {B} to an empty box, by setting all coordinate intervals {C[0..d-1]}
    to some empty interval. */

void box_include_point(box_dim_t d, interval_t B[], double p[], interval_t C[]);
  /* Sets box {C[0..d-1]} to the smallest box enclosing box {B} and point {p[0..d-1]}.
    In particular, if {B} is empty, the result is the trivial box that has only {p}. */

void box_join(box_dim_t d, interval_t A[], interval_t B[], interval_t C[]);
  /* Sets box {C[0..d-1]} to the smallest box enclosing both boxes {A} and {B}.]
    If either box is empty, the result is the other box. */

void box_meet(box_dim_t d, interval_t A[], interval_t B[], interval_t C[]);
  /* Sets box {C[0..d-1]} to the intersection of boxes {A} and {B}.
    Note that the result may be empty even if {A} and {B} are not. */

void box_widen(box_dim_t d, interval_t B[], double margin, interval_t C[]);
  /* Sets box {C[0..d-1]} to the box {B} with the limits on each axis
    displaced outwards by {margin}. The displacement will be inwards if
    {margin} is negative; in this case the box {C} may be empty. If {B}
    is empty, {C} will be empty. */

void box_widen_var(box_dim_t d, interval_t B[], double mrglo[], double mrghi[], interval_t C[]);
  /* Sets box {C[0..d-1]} to the box {B} with the low limit on each axis
    {i} decremended by {mrglo[i]}, and the high limit on each axis {i}
    incremented by {mrghi[i]}. If any of these displacements is
    negative, the box may become empty. */

void box_round(box_dim_t d, interval_t B[], double unit, bool_t out, interval_t C[]);
  /* Sets box {C[0..d-1]} to the box {B} with the limits along each axis
    {i} rounded to the nearest multiple of {step}. The limits are
    rounded outwards if {out} is true, inwards if {out} is false. */

void box_split
  ( box_dim_t d, 
    interval_t B[],  
    box_axis_t a,  
    double x,  
    interval_t BLO[],  
    interval_t BMD[],
    interval_t BHI[] 
  );
  /* Splits the box {B[0..d-1]} by a plane perpendicular to axis {a}
    (which must be in {0..d-1}) at coordinate {x}, into boxes
    {BLO[0..d-1]} (below {x}), {BMD[0..d-1]} (at {x}) and {BHI[0..d-1]}
    (above {x}).
    
    If {B} is empty, all three boxes {BLO,BMD,BHI} will be empty. 
    
    Otherwise, two of the boxes may be empty, depending on the position
    of {x} relative to the interval {B[a]} and on whether {B[a]} is a
    singleton or an open interval.
    Specifically,
    
      {BLO} will be empty iff {x <= LO(B[a])};
      
      {BHI} will be empty iff {x >= HI(B[a])};
      
      {BMD} will not be empty iff {LO(B[a] < x < H(B[a])}
            or {LO(B[a]) == x == HI(B[a])}
    
    if {LO(B[a]) < x < HI(B[a]}, the three boxes will be not empty,
    and will have {BLO[a] = (LO(B[a]) _ x)}, {BMD[a] = [x_x]}, and {BHI = (x, BHI(a))}.
    Otherwise the only non-empty box will be the whole of {B}. In any case, the result boxes that
    are not empty will be equal to {B} along all axes except {a}.
    
    If any of the parameters {BLO,BMD,BHI} is {NULL}, the corresponding
    box is not generated. */

void box_throw(box_dim_t d, double elo, double ehi, double p_empty, double p_single, interval_t B[]);
  /* Stores into {B[0..d-1]} a random {d}-dimensional box. 
  
    With probability {p_empty}, the box will be empty.
    
    If the box is not empty, every projection {B[i]} will be strictly
    contained in the open interval {(elo _ ehi)}, which must not be
    empty. Each {B[i]} will be a singleton in that range, independently,
    with probability {p_single}. Otherwise it will be an non-empty,
    non-trivial sub-intrerval of that range. */

/*
  BOX MAPPING AND UNMAPPING
  
  A non-empty {m}-dimensional box {B[0..d-1]} with spanning axes {s[0..m-1]}
  and stabbing axes {t[0..d-m-1]} defines a mapping from {R^m} to
  {R^d}, whereby each coordinate {z[j]} of the argument is scaled and
  shifted so that the interval [0_1] becomes the interval {x[s[j]] =
  B[s[j]]}, and {x[t[j]]} is set to the singleton {B[t[j]}.
  
  We will denote this transformation by {x = S_B(z)}. The numbers
  {z[0..m-1]} are the /coordinates of {x} relative to {B}/. */

void box_point_map(box_dim_t d, double z[], interval_t B[], double x[]);
  /* Computes the point of of {R^d} which has coordinates {z[0..m-1]}
    relative to the {m}-dimensional box {B[0..d-1]}, and stores it in
    {x[0..d-1]}. Fails if {B} is empty. */

void box_point_unmap(box_dim_t d, double x[], interval_t B[], double z[]);
  /* Computes the cordinates of point {x[0..d-1]} relative to the
    {m}-dimensional box {B}, and stores them in {z[0..m-1]}. The point
    {x} is implicitly projected onto the affine span of {B}. */

void box_box_map(box_dim_t d, interval_t Z[], interval_t B[], interval_t X[]);
  /* Computes the box of {R^d} which has coordinate ranges {Z[0..m-1]}
    relative to the {m}-dimensional box {B[0..d-1]}, and stores it
    in {X[0..d-1]}. */

void box_box_unmap(box_dim_t d, interval_t X[], interval_t B[], interval_t Z[]);
  /* Computes the coordinate ranges of box {X[0..d-1]} relative to the
    {m}-dimensional box {B[0..d-1]}, and stores them in {Z[0..m-1]}.
    The box {X} is implicitly projected onto the affine span of {B}. */

/* 
  BOX ORIENTATION
  
  A box {B} is/parallel/ to a coordinate axis {i} if its projection on
  that axis is a non-trivial interval; otherwise {B} is
  /perpendicular/ to axis {i}. We also say that axis {i} is /spanning/
  or /normal/ to the box, respectively.
  
  We denote by {Spn(B)} the spanning axes of box {B}, and by {Nrm(B)}
  the normal ones. The set {Nrm(B)} is the /orientation/ of the
  box. */
  
typedef set32_t box_axis_set_t;
  /* A set of axes, e.g. spanning or normal axes of a box. */
  
typedef set32_index_t box_axis_index_t;
  /* The index of an axis in a set of axes. */

#define box_NO_AXIS_INDEX (set32_NO_INDEX)

/* 
  FACETS AND FACES
  
  A /facet/ of an {m}-dimensional box {B} is a maximal box with
  dimension {m-1} contained in the boundary of {B}. There are {2*m}
  facets, two for each axis parallel to {B}: the /inferior/ facet, in
  the {-oo} direction from {B}, and the /superior/ one, in the {+oo}
  direction.
  
  The /faces/ of {B} are the box {B} itself, and the facets of all its
  faces. The /inferior/ (resp. /superior/) faces are the inferior
  (resp. superior) facets of all faces of {B}. An {m}-dimensional box
  has {3^m} faces, of which {2^d-1} are inferior and {2^d-1} are
  superior.
  
  RELATIVE FACE POSITIONS
  
  An {m}-dimensional box {B} has {3^m} faces. Let {Spn(B) = ax[0..m-1]}
  be the spanning axes of {B}. Each face {F} of {B} can be identified by its
  /signature/, a vector {dir[0..m-1]} where {dir[j]} specifies the
  projection {F_j} of {F} onto axis {ax[j]}, in terms of the
  projection {B_j} of {B} on that axis. Namely,
  
    {F_j = lo(B_j)}  if {dir[j] = -1}
    {F_j = B_j}      if {dir[j] = 0}, and, 
    {F_j = hi(B_j)}  if {dir[j] = +1}.
  
  In particular, if {dir = (0,.. 0)}, then {F} is {B} itself. An
  inferior face of a box has a signature consisting entirely of {0}s
  and {-1}s; the corresponding superior face has the {-1}s replaced by
  {+1}s. The signatures {(-1,.. -1)} and {(+1,.. +1)} specify the
  inferior and superior corners of {B}; respectively.
  
  */
 
typedef enum { box_NEG = -1, box_ZER = 0, box_POS = +1 } box_signed_dir_t;
  /* Three-valued direction along an axis: {NEG} towards {-oo}, {POS}
    towards {+oo}, {ZER} in the middle/interior. Used to specify the
    position of a face relative to a box, etc. */
  
/* 
  RELATIVE FACE INDICES
  
  The signature {dir[0..m-1]} of a face {F} relative to a box {B}
  can be packed into a single /relative face index/ {fi}, by the
  formula {fi = SUM{ 3^j * (dir[j] \bmod 3) : j = 0..m-1}}.
  
  Note that that the relative face index {fi} ranges from 0 to {3^m-1}.
  The face index {fi = 0} identifies the box {B} itself.
  The inferior corner has index {fi = (3^{m}-1)/2}; and the
  superior corner has index {fi = 3^{m}-1}. */
  
typedef uint32_t box_face_index_t;
  /* Identifies one of the {3^m} faces of an {m}-dimensional box. */

box_dim_t box_face_dimension(box_dim_t m, box_face_index_t fi);
  /* Returns the dimension the face of an {m}-dimensional 
    box {B} that is identified by the relative face index {fi}. */

box_signed_dir_t box_face_position(box_face_index_t fi, box_axis_index_t j);
  /* Given the index {fi} of a face {F} of some box {B}, returns
    the position of {F} relative to {B} along axis {Spn(B)[j]}. */

void box_face_signature(box_dim_t m, box_face_index_t fi, box_signed_dir_t dir[]);
  /* Given the index {fi} of a face {F} of an {m}-dimensional
    box {B}, stores in {dir[j]} the position of {F} relative
    to {B} along axis {Spn(B)[j]}. */

void box_face_rel_box(box_dim_t m, box_face_index_t fi, interval_t F[]);
  /* Given the relative index {fi} of an {m}-dimensional face {F} of a {d}-dimensional
    box {B}, stores in each {F[j]} the coordinate range of face
    {fi} along axis {Spn(B)[j]}, for {j = 0..m-1}.
    
    The coordinates are relative to the box {B}, i.e. 0 means the inferior
    side of the box, 1 means the superior side. Thus, {F[j]} will be
    set to either {0,0}, {1,1}, or {0,1}. Note that the first two mean
    singleton sets {\set{0}} and {\set{1}}, while the latter means the
    open interval {(0_1)}. */

void box_print(FILE *wr, box_dim_t d, interval_t B[]);
  /* Prints the box {B} as a series of intervals, in some default format */
  
void box_gen_print
  ( FILE *wr, 
    box_dim_t d, 
    interval_t B[], 
    char *fmt, 
    char *pref, 
    char *sep, 
    char *suff
  );
  /* Prints the box {B} as a series of intervals, preceded by {pref},
    separated by {sep}, and followed by {suff}. Each end of the interval
    is printed with format {fmt}, which must have exactly one "%f",
    "%g", or "%e" spec. 
    
    If the box is empty, all intervals are printed as "[]". If an
    interval is a single value ({LO == HI}), it is printed as that value
    between "{" and "}". Otherwise each interval is printed as
    "({LO}_{HI})" where {LO} and {HI} are the endpoints.
    
    Suitable defaults are provided if any of {fmt,pref,sep,suff} are
    {NULL}. */

void box_face_print(FILE *wr, box_dim_t m, box_face_index_t fi);
  /* Prints the signature of face number {fi} of an {m}-dimensional box. */
  
#endif
