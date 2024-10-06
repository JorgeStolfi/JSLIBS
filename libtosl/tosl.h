/* Basic defs for the topological slicing algorithm of Minetto et al. (2024).  */
/* Last edited on 2024-10-06 16:46:37 by stolfi  */

#ifndef tosl_H
#define tosl_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

typedef int32_t tosl_arc_id_t;
  /* Each arc (oriented edge) {a} of the mesh is represented in this
    algorithm by an integer index {ia} into the arc table
    {Arc[0..2*NE-1]}. A {tosl_arc_id_t} value is such an index.

    In terms of {haf.h}, the index {ia} is {haf_arc_id(a)}. The value
    {-1} is used to signify "no such arc" (that is, {a} is {NULL}). */

typedef int32_t tosl_vert_id_t;
  /* Each vertex {v} of the mesh is represented in this algorithm by an
    integer index {iv} into the table {Vpos[0..NV-1]}. A
    {tosl_vert_id_t} value is such an index. */

typedef int32_t tosl_plane_id_t;
  /* Type of index of a slicing plane. */

typedef int32_t tosl_coord_t;
  /* Quantized coordinate (integer multiple of the fund. unit {\lambda}. */

typedef struct tosl_point_t { tosl_coord_t c[3]; } tosl_point_t;
  /* Quantized vertex coords. */

typedef struct tosl_arc_t 
  { tosl_arc_id_t skip;   /* The {.lnext} link, with shortcuts. */
    tosl_vert_id_t ivorg; /* Index of origin vertex. */
    /* Links for active/bucket arc lists:  */
    tosl_arc_id_t pred;   /* Index of previous edge rec in arc list. */
    tosl_arc_id_t succ;   /* Index of next edge rec in arc list. */
  } tosl_arc_t;
  /* A {tosl_arc_t} record contains the data about an oriented arc
    of the mesh. It corresponds to an {haf_arc_t} of {haf.h}

    The field {Arc[ia].skip} represents the auxiliary link {a.skip} of
    the topological slicing algorithm as described in the paper. See
    {topolic_mesh_slice.h} below.

    For each arc {a} with id {ia}, {Vpos[Arc[ia].ivorg]} are the
    coordinates of the origin vertex of {a}.

    The fields {.pred} and {.succ} are indices into {Arc} that are used
    to link the {tosl_arc_t} records into circular doubly-linked
    lists that comprise the arc bucket lists an the active arc lists of
    the algorithm. If {Arc[ia]} is not part of any such list, then
    {Arc[ia].pred} and {Arc[ia].succ} are both {-1}. */

#define tosl_sym(ia) ((ia) ^ 1)
  /* Given the index {ia} of an arc {a}, returns the index of the arc {haf_sym(a)}. */
 
/* MISCELLANEOUS */

char *tosl_arc_id_to_string(tosl_arc_id_t ka);
  /* Converts the arc identifier {ka} to a string "a{e}:{o}" where
    {e} is the edge id (namely {ka/2}) and {o} is the orientation
    bit (namely {ka%2}).  However, if {ka} is {-1}, returns "*:*"
    instead.  The result is always a newly allocated string. */
        
void tosl_arc_id_print(FILE *wr, char *pref, tosl_arc_id_t ka, char *suff);
  /* Prints to {wr} the arc id {ka} formatted by {tosl_arc_id_to_string}, preceded 
    by {pref} and followed by {suff}. */
        
void tosl_tri_arc_id_print(FILE *wr, char *pref, tosl_arc_id_t ka0, tosl_arc_id_t ka1, tosl_arc_id_t ka2, char *suff);
  /* Prints to {wr} the arc ids {ka0,ka1,ka2}, each formatted by {tosl_arc_id_to_string},
    separated by " - ", preceded by {pref} and followed by {suff}. */
  
double tosl_user_cpu_time_usec(void);
  /* Current user-mode CPU time accumulated by process. */

void tosl_compute_avg_dev(int32_t NT, double time[], double *avg_P, double *dev_P);
  /* Sets {*avg_P} and {*dev_P} to the average and deviation of {time[0..NT-1]} */
  
#endif
