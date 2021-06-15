/* r3_motion.h --- paths in R^3 with local frames */
/* Last edited on 2021-06-12 09:27:27 by jstolfi */

#ifndef r3_motion_H
#define r3_motion_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3x3.h>
#include <vec.h>

/* PATH STATES */

typedef struct r3_motion_state_t 
  { r3x3_t M;  /* Scaling and orientation matrix. */
    r3_t p;    /* Position of reference point. */
  } r3_motion_state_t;
  /* The position, scaling, and orientation of an object in space.
    
    Each point {q} of the original reference object is assumed to be
    mapped by {r3x3_map_row(q,.M)} and then translated by adding {.p}.
    
    The rows of the matrix {.M} can be seen as three vectors of a local
    coordinate frame at {S.p}. They are referred to as {S.u}, {S.v}, and {S.w}.
    
    For an Euclidean transformation, {.M} should be orthogonal, with all
    three row vectors pairwise orthogonal and with the same length. For
    an isometry, in particular, {.M} should be orthonormal, i.e. the
    vectors {S.u}, {S.v}, and {S.w} shoudl have unit length. */
    
typedef r3_motion_state_t r3_motion_proc_t(double t);
  /* Type of a function that describes the motion and deformation of
    an object in 3-space. Both the {M} and {p} components of the state must
    be continuous functions of {t}. */

void r3_motion_state_canonical(r3_motion_state_t *C);
  /* Stores into {C} the canonical state, whose position is the origin
    and whose pose matrix is the identity. */

void r3_motion_state_compose(r3_motion_state_t *S, r3_motion_state_t *T, r3_motion_state_t *C);
  /* Stores in state {C} a state that is related to {S} as {T} is related to the 
    origin.  
    
    Namely, there is a unique isometry of {\RR^3} that takes the 
    canonical state to state {T}; namely, rotate by {T.M} about the
    origin, then translate by {T.p}.  The procedure applies that same 
    isometry to {S}, and stores the result in {C}. */

/* PATH SAMPLING */

void r3_motion_sample_uniform
  ( r3_motion_proc_t path, 
    double t0,
    double t1,
    int32_t n,
    bool_t mids,
    double t[],
    r3_motion_state_t S[]
  );
  /* Samples {path(t)} for {n} equally spaced values of {t} between
    {t0} and {t1}.  Returns those values of {t} and the corresponding
    path states in the arrays {t[0..n-1]} and {S[0..n-1]}, if they are
    not {NULL}. 
    
    If {mids} is false, will set {t[0]=t0} and {t[n-1]=t1}. In that
    case, {n} must be at least 2, and the sampling step {dt} will be
    {(t1-t0)/(n-1)}.
    
    If {mids} is true, the step {dt} will be {(t1-t0)/n}, and the
    procedure will set {t[0]=t0+dt/2}, {t[n-1]=t1-dt/2}. This may be the
    appropriate option when the path is closed (i.e.
    {path(t0)==path(t1)}). */

/* STANDARD PATHS

  Each function in this section computes a state {S=S(t)} along some
  standard path. The orientation matrix {S.M} is an orthonormal frame. If the
  curve is has a well-defined velocity vector, then {S.u} (row 0 of {S.M}) is its
  direction, a unit vector tangent to the curve. If the curve is not
  straight, {S.v} (row 1 of {S.M}) points directly towards the center of
  curvature; {S.u} and {S.v} then span the osculating plane of the curve at
  that point. Then {S.w} (row 2 of {S.M}) is the normal of that plane, oriented
  by the right-hand rule. If the curve has zero or discontinuous
  velocity, or is straight, the frame is defined conventionally for each
  case. */

void r3_motion_circle(double t, double L, double A, r3_motion_state_t *S);
  /* Stores in {S} a state at time {t} along a circle
    on the coordinate plane {Z=0}. 
    
    The point {S.p} starts at the origin when {t=0}; when {t=1}, {S.p}
    has traversed a distance {L} along the circle and turned by an angle
    {A} (in radians, counterclockwise). If {L>0}, the arc starts tangent
    to the {X} coordinate axis in the direction {+X}.  The parameter
    {t} can extend beyond 
    
    The parameters {L} and {A} can be zero or negative. In particular,
    if {A=0} the arc is a straight line that ends at {(L,0,0)}. If {L<0}
    the arc starts moving in the {-X} direction; however, the pose when
    {t=0} does not depend on the signs of {L} or {A}. */

void r3_motion_helix(double t, double L, double A, double H, r3_motion_state_t *S);
  /* Stores in {S} a state at time {t} along a helix curve. 
  
    The helix lies on a cylindrical surface whose axis is parallel to
    the {Z} coordinate axis. The point {S.p} starts at the origin when
    {t=0}. The projection of the helix on the coordinate plane {Z=0} is
    the circle described by {r3_motion_circle} with the same parameters
    {L} and {A}. The difference is that the {Z} coordinate of the point
    {S.p} increases at an uniform rate, and is {H} when {t=1}.
    
    In particular, if {H} is zero, the result is identical to that of
    {r3_motion_circle(t,L,A)}. If {H} is not zero, rows 0 and 2 of {S.M}
    are rotated around row 1 so that row 0 is tangent to the curve. */

void r3_motion_state_print (FILE *f, r3_motion_state_t *S);
  /* Prints {S} on file {f}, with some default format. */

void r3_motion_state_gen_print 
 ( FILE *f,
   r3_motion_state_t *S, 
   char *pfmt, 
   char *Mfmt, 
   char *olp, char *osep, char *orp, /* Outer delimiters for matrix. */
   char *ilp, char *isep, char *irp  /* Delimiters for matrix rows and position. */
  );
  /* Prints position {S.p} and the three rows of the orientation matrix
    {S.M} on file {f}, formatting each coordinate with {pfmt} and
    {Mfmt}, respectively. 
    
    The strings {olp} and {orp} are printed respectively before the
    first and after the last of these four vectors; the string {osep} is
    printed between them. The letters "p", "u", "v", "w" are printed
    before each vector {x}, and then the vector is printed using
    {r3_gen_print(f,x,fmt,ilp,isep,irp)}. Defaults are provided for any
    of these strings which are NULL. */
 
void r3_motion_state_debug(FILE *wr, r3_motion_state_t *S, char *indent, char *title);
  /* Prints state {S} for debugging, with the given identation and title strings. */
   
#endif
