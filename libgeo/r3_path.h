/* r3_path.h --- Paths defined by points, velocities, times */
/* Last edited on 2024-12-05 10:28:02 by stolfi */

#ifndef r3_path_H
#define r3_path_H

#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>

/* !!! Move {r3_path_state_t} to {r3_motion}, rename {r3_motion} to {r3_motion}. !!! */

typedef struct r3_path_state_t 
  { double t;  /* Nominal time. */
    r3_t p;    /* Position. */
    r3_t v;    /* Velocity. */
  } r3_path_state_t;
  /* Parameters of state along a path, consisting of a time value {t},
    a point {p}, and a velocity vector {v}. */

r3_path_state_t r3_path_state_from_r3_motion_state(double t, r3_motion_state_t *S, double v);
   /* Creates a {r3_path_state_t} {T} with time {t}, position {S.p} and moving with speed {v}
     in the direction of of {S.u} (row 0 of {S.M}). */ 

r3_motion_state_t r3_path_state_to_r3_motion_state(r3_path_state_t *S);
   /* Creates n {r3_motion_state_t} {T} from the given {r3_path_state_t} {S}. 
     The position {T.p} is {S.p}, and the first vector {T.u} ({=T.M[0][*]}) is the 
     unit direction vector of the velocity vector {S.v}.  The 
     second vector {T.v} ({=T.M[1][*]}) is an arbitrary vector orthogonal to
     {T.u}, and the third vector {T.w} ({=T.M[2][*]}) is the cross product of the two. */ 

r3_path_state_t r3_path_state_from_r3_motion_state(double t, r3_motion_state_t *S, double v);
   /* Creates a {r3_path_state_t} with time {t}, position {S.p} and moving with speed {v}
     in the direction of the {Z} axis of {S} (row 2 of {S.M}. */ 

double r3_path_length_estimate(r3_path_state_t *S, r3_path_state_t *T, int32_t order);
  /* Estimates the length of the path from state {S} to state {T}
     by converting it to a Bezier arc with {r3_path_bezier_from_states}
     and calling {r3_bezier_length_estimate}. */

void r3_path_interpolate_some(r3_path_state_t S[], uint32_t N);
  /* Computes the times, positions and velocities of {S[1..N-2]} to
    smoothly interpolate between {S[0]} and {S[N-1]},
    which are not changed. Currently uses the {r3_path_bezier_from_states}
    and {r3_bezier_split}. */
   
void r3_path_bezier_from_states(r3_path_state_t *S, r3_path_state_t *T, r3_t *p1, r3_t *p2);
  /* Computes the two middle Bezier control points of a cubic arc that starts 
    at time {S.t} and position {S.p} with velocity {S.v} and ends at time {T.t} 
    and position {T.p} with velocity {T.v}. The control points
    of the arc will be {p0=S.p}, {p1}, {p2}, and {p3=T.p}. */
 
/* PRINTOUT */

void r3_path_state_print (FILE *f, r3_path_state_t *S);
  /* Prints {S} on file {f}, with some default format. */

void r3_path_state_gen_print
 ( FILE *f,
   r3_path_state_t *S,
   char *tfmt,
   char *pfmt,
   char *vfmt,
   char *olp, char *osep, char *orp, /* Outer delimiters for matrix. */
   char *ilp, char *isep, char *irp  /* Delimiters for matrix rows and position. */
  );
  /* Prints the time {S.t}, the position {S.p}, and the velocity {S.v}
    on file {f}, formatting each coordinate with {tfmt}, {pfmt}, and
    {vfmt}, respectively.
    
    The strings {olp} and {orp} are printed respectively before {S.t}
    and after {S.v}; the string {osep} is printed between the three
    components. The letters "t", "p", and "v" are printed before each
    component. The time is then printed with format {tfmt} bracketed by
    {ilp} and {irp}. Each vector component {x} then the vector is
    printed using {r3_gen_print(f,x,fmt,ilp,isep,irp)}. Defaults are
    provided for any of these strings which are NULL. */
  
void r3_path_state_debug(FILE *wr, r3_path_state_t *S, char *indent, char *title);
  /* Prints state {S} for debugging, with the given identation and title strings. */

#endif

