/* voxm_path.h --- path operations for voxel-based modeling */
/* Last edited on 2016-04-04 00:10:54 by stolfilocal */

#ifndef voxm_path_H
#define voxm_path_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_path.h>

/* !!! Move {voxm_path_state_t} to {r3_path}, rename {r3_path} to {r3_motion}. !!! */

typedef struct voxm_path_state_t 
  { double t;  /* Nominal time. */
    r3_t p;    /* Position of fork. */
    r3_t v;    /* Velocity. */
  } voxm_path_state_t;
  /* Parameters of state along a path, consisting of a time value {t},
    a point {p}, and a velocity vector {v}. */

voxm_path_state_t voxm_path_state_from_r3_path_state(double t, r3_path_state_t *S, double v);
   /* Creates a {voxm_path_state_t} {T} with time {t}, position {S.p} and moving with speed {v}
     in the direction of of {S.u} (row 0 of {S.M}). */ 

r3_path_state_t voxm_path_state_to_r3_path_state(voxm_path_state_t *S);
   /* Creates n {r3_path_state_t} {T} from the given {voxm_path_state_t} {S}. 
     The position {T.p} is {S.p}, and the first vector {T.u} ({=T.M[0][*]}) is the 
     unit direction vector of the velocity vector {S.v}.  The 
     second vector {T.v} ({=T.M[1][*]}) is an arbitrary vector orthogonal to
     {T.u}, and the third vector {T.w} ({=T.M[2][*]}) is the cross product of the two. */ 

voxm_path_state_t voxm_path_state_from_r3_path_state(double t, r3_path_state_t *S, double v);
   /* Creates a {voxm_path_state_t} with time {t}, position {S.p} and moving with speed {v}
     in the direction of the {Z} axis of {S} (row 2 of {S.M}. */ 

double voxm_path_length_estimate(voxm_path_state_t *S, voxm_path_state_t *T, int order);
  /* Estimates the length of the path from state {S} to state {T}
     by converting it to a Bezier arc with {voxm_bezier_from_path_states}
     and calling {voxm_bezier_length_estimate}. */

void voxm_path_interpolate_some(voxm_path_state_t S[], int N);
  /* Computes the times, positions and velocities of {S[1..N-2]} to
    smoothly interpolate between {S[0]} and {S[N-1]},
    which are not changed. Currently uses the {voxm_bezier_from_path_states}
    and {voxm_bezier_split}. */
    
/* PRINTOUT */

void voxm_path_state_print (FILE *f, voxm_path_state_t *S);
  /* Prints {S} on file {f}, with some default format. */

void voxm_path_state_gen_print
 ( FILE *f,
   voxm_path_state_t *S,
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
  
void voxm_path_state_debug(FILE *wr, voxm_path_state_t *S, char *indent, char *title);
  /* Prints state {S} for debugging, with the given identation and title strings. */

#endif

