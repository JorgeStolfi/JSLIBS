/* See {tosl_mesh_slice.h} */
/* Last edited on 2024-10-09 10:18:17 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>

#include <haf.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_arc_list.h>
#include <tosl_slice.h>

#include <tosl_mesh_slice.h>

void tosl_extract_contour
  ( tosl_arc_id_t ia,
    tosl_slice_t *S,
    tosl_arc_id_t *A_P,
    tosl_arc_id_t *B_P,
    tosl_mesh_t *mesh,
    int8_t debug
  );
  /* Given the index {ia} of an upward arc {a} that crosses the current slicing plane,
    repeatedly avances {ia} to the index of the first (and only) upward edge after {haf_sym(a)}, in 
    CCW in the right face of {a}, that crosses
    the slicing plane {S.Zp}; until coming back to the starting arc {a}.  Appends the arcs
    found in that traversal to {S.iarc}, marking the last one as end of another contour. 
    Also removes those arcs from the active list {*A_P} and moves them to the auxiliary
    list {*B_P}. */

void tosl_mesh_slice
  ( tosl_mesh_t *mesh,
    int32_t NP,
    tosl_coord_t Zplane[],
    tosl_arc_id_t L[],
    tosl_post_proc_t *post,
    double *time_P,
    int8_t debug
  )
  {
    if (debug) { fprintf(stderr, "> --- %s --------------------\n", __FUNCTION__); }

    tosl_arc_id_t B = -1; /* Active list remaining from previous iteration. */
    
    if (time_P != NULL) { (*time_P) = 0.0; }
    
    for (tosl_plane_id_t ip = 0; ip < NP; ip++)
      { 
        if (debug) { fprintf(stderr, "  plane %d Z = %d", ip, Zplane[ip]); }
        double tstart, tstop;
        if (time_P != NULL) { tstart = tosl_user_cpu_time_usec(); }
        
        tosl_coord_t Zp = Zplane[ip];
        /* Get the active list {A} for this iteration. */
        tosl_arc_id_t A = tosl_arc_list_merge(&B, &(L[ip]), Zp, mesh);
        assert(B == -1); 
        
        /* Allocate the slice storage area: */
        int32_t NA = tosl_arc_list_len(A, mesh); 
        tosl_slice_t *S = tosl_slice_new(NA, Zp);
        if (debug) { fprintf(stderr, " |A| = %d\n", NA); }
        
        while (A != -1)
          { /* Pop an active edge: */
            tosl_arc_id_t ia = A;
            tosl_extract_contour(ia, S, &A, &B, mesh, debug);
          }
        
        if (time_P != NULL) { tstop = tosl_user_cpu_time_usec();  (*time_P) += (tstop - tstart); }
        
        post(ip, S);
      }
    if (debug) { fprintf(stderr, "< --- %s --------------------\n", __FUNCTION__); }
  }  

void tosl_extract_contour
  ( tosl_arc_id_t ia,
    tosl_slice_t *S,
    tosl_arc_id_t *A_P,
    tosl_arc_id_t *B_P,
    tosl_mesh_t *mesh,
    int8_t debug
  )
  {
    tosl_coord_t Zp = S->Z;

    if (debug) { fprintf(stderr, "  > --- %s --------------------\n", __FUNCTION__); }
    if (debug) { fprintf(stderr, "    Zp = %+d\n\n", Zp); }
    if (debug) { tosl_mesh_arc_print(stderr, "    ia = ", ia, "\n\n", mesh); }
    
    tosl_arc_id_t ka = ia;
    do
      {
        if (debug) { tosl_mesh_arc_print(stderr, "      ka = ", ka, "\n", mesh); }
    
        /* Get {Z} of origin and destination of {ka}:  */
        tosl_coord_t Zorg = mesh->Vpos[mesh->Arc[ka].ivorg].c[2]; /* {Z} of origin of {ka}.  */
        tosl_arc_id_t sa = tosl_sym(ka);
        tosl_coord_t Zdst = mesh->Vpos[mesh->Arc[sa].ivorg].c[2]; /* {Z} of destination of {ka}.  */

        /* At this point {Arc[ka]} should be an upward-pointing arc that crosses the plane.  */
        assert(Zorg < Zp);
        assert(Zdst > Zp);

        /* Append to the slice the intersection of {ka} with the plane:  */
        assert(S->NV < S->NV_max);
        S->iarc[S->NV] = ka;
        (S->NV)++;

        /*  Move arc {ka} from set {A} to set {B} */
        tosl_arc_list_remove(A_P, ka, mesh);
        tosl_arc_list_add(B_P, ka, mesh);

        /* Search for next up-crossing arc {ja} using skip pointers:  */
        tosl_arc_id_t ja = sa;
        while (1) { 
          ja = mesh->Arc[ja].skip;
          if (debug) { tosl_mesh_arc_print(stderr, "        ja = ", ja, "\n", mesh); }
          tosl_vert_id_t jv = mesh->Arc[tosl_sym(ja)].ivorg;
          if (mesh->Vpos[jv].c[2] > Zp) { break; }
        }
        
        /* Update the skip pointer to short-cut from {ka} directly to {ja}:  */
        if (debug) 
          { tosl_arc_id_print(stderr, "      updating ", sa, ".skip");
            tosl_arc_id_print(stderr, " to ", ja, "\n");
          }
        mesh->Arc[sa].skip = ja;
        
        /* Prepare for next iteration:  */
        ka = ja;
      }
    while (ka != ia);
    /* Mark the {S->iarc[S->NV-1]} as the last vertex of a contour:   */
    S->iarc[S->NV-1] = tosl_sym(S->iarc[S->NV-1]);

    if (debug) { fprintf(stderr, "  < --- %s --------------------\n", __FUNCTION__); }
    return;
  }
