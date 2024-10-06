/* See {tosl_build_lists_bin_sec.h}. */
/* Last edited on 2024-10-06 16:48:41 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_arc_list.h>
#include <tosl_build_lists_bin_sec.h>

tosl_plane_id_t tosl_build_lists_bin_sec_search
  ( tosl_coord_t Ze,
    int32_t NP,
    tosl_coord_t Zplane[],
    int32_t use_bin,
    int32_t use_sec,
    int32_t debug,
    int32_t *ns_P
  );
  /* Given: the lower Z coordinate {Ze} of an edge, and the list of plane
    Z-coordinates {Zplane[0..NP-1]}. If {Ze < Zplane[0]}, returns 0;
    if {Ze > Zplane[NP-1], returns {NP}; otherwise returns an index
    {i} in {1..NP-1} such that {Zplane[i-1] < Ze < Zplane[i]}.
    Assumes that the array Zplane[] is sorted in strictly increasing
    order. Assumes that {Ze} is an even integer and all {Zplane[i]}
    are odd integers.

    The parametes {use_bin} and {use_sec} must be 0 (false) or 1 (true).
    If {use_bin} is true, uses a binary search step at each iteration.
    If {use_sec} is true, uses a secant step at each iteration.
    If both are true, alternates between the two steps. 
    They must not be both false.
    
    If {ns_P} is not {NULL}, adds to {*ns_P} the number of binary and/or secant steps
    performed.  Note that the caller must initialize {*ns_P} before the first call. */

tosl_arc_id_t *tosl_build_lists_bin_sec
  ( int32_t NP,
    tosl_coord_t Zplane[],
    tosl_mesh_t *mesh,
    int32_t use_bin,
    int32_t use_sec,
    int32_t debug
  )
  {
    assert((mesh->NE > 0) && (mesh->NE <= tosl_mesh_MAX_EDGES));
    assert((NP > 0) && (NP <= tosl_mesh_MAX_PLANES));
 
    int32_t NA = 2*mesh->NE;
    
    /* Allocate the bucket table and set all entries to empty:  */
    tosl_arc_id_t *L = malloc(NP*sizeof(tosl_arc_id_t));
    assert(L != NULL);
    for (tosl_plane_id_t ip = 0; ip < NP; ip++) { L[ip] = -1; }
    
    int32_t ns = 0;
    int32_t *ns_P = (debug ? &ns : NULL);
    
    for (tosl_arc_id_t ia = 0; ia < NA; ia++)
      { /* Make {Arc[ia]} into a singleton list: */
        mesh->Arc[ia].pred = ia;
        mesh->Arc[ia].succ = ia;
        
        /* Get the origin and destination {Z}-coordinates {Zorg,Zdst}: */
        tosl_coord_t Zorg = mesh->Vpos[mesh->Arc[ia].ivorg].c[2]; /* {Z} of origin of {ia}.  */
        tosl_arc_id_t ja = tosl_sym(ia);
        tosl_coord_t Zdst = mesh->Vpos[mesh->Arc[ja].ivorg].c[2]; /* {Z} of destination of {ia}.  */
        
        /* Add the arc to a list {L[ip]} if appropriate: */
        if ((Zorg < Zdst) && (Zorg < Zplane[NP-1]) && (Zdst > Zplane[0]))
          { /* Arc is up-going and is not above or below all planes. */
            tosl_plane_id_t ip = tosl_build_lists_bin_sec_search
              ( Zorg, NP, Zplane, use_bin, use_sec, debug, ns_P );
            /* Paranoia, checking {tosl_plane_search}: */
            assert(Zorg < Zplane[ip]);
            assert((ip == 0) || (Zorg > Zplane[ip-1]));
            if (Zdst > Zplane[ip]) 
              { /* Arc {ia} crosses at least one slicing plane: */
                tosl_arc_list_add(&(L[ip]), ia, mesh);
              }
            else
              { /* Check that no vertices are on planes: */
                assert(Zdst < Zplane[ip]); 
              }
          }
        else
          { /* Check that no vertices are on planes: */
            assert((Zorg >= Zdst) || (Zorg > Zplane[NP-1]) || (Zdst < Zplane[0]));
          }
      }
    if (debug)
      { assert(ns_P != NULL);
        double nspe = ((double)ns)/((double)mesh->NE);
        fprintf(stderr, "performed %.2f binary/secant steps per edge\n", nspe);
      }
    return L;
  }

tosl_plane_id_t tosl_build_lists_bin_sec_search
  ( tosl_coord_t Ze,
    int32_t NP,
    tosl_coord_t Zplane[],
    int32_t use_bin,
    int32_t use_sec,
    int32_t debug,
    int32_t *ns_P
  )
  {
    if (debug) { fprintf(stderr, "  Ze = %+10d, bin:%c sec:%c\n", Ze, "FT"[use_bin], "FT"[use_sec]); } 
    if (Ze < Zplane[0]) { return 0; }
    if (Ze > Zplane[NP-1]) { return NP; }
    assert(use_bin || use_sec);
    tosl_plane_id_t imin = 0, imax = NP-1;
    while (1)
      { /* At this point {imin < imax} and {Zplane[imin] < Ze < Zplane[imax]}. */
        if (imax - imin <= 1) { break; }
        if (use_sec && (imax - imin >= 3))
          { /* Do a secant search step: */
            int64_t dzedge = Ze - Zplane[imin];
            int64_t dzplane = Zplane[imax] - Zplane[imin];
            tosl_plane_id_t isec = imin + (int32_t)((dzedge*(imax - imin))/dzplane);
            /* Pull {isec} towards midpoint:  */
            tosl_plane_id_t ictr = (imin + imax)/2;
            int32_t di = (isec - ictr)/64;
            isec = isec - di;
            /* Keep {isec} away from {imin,imax}: */
            if (isec <= imin) { isec = imin + 1; }
            if (isec >= imax) { isec = imax - 1; }
            assert((imin < isec) && (isec < imax));
            /* Split the range at {isec}:  */
            if (debug) { fprintf(stderr, "  %9d .. %9d  %5.2f sec: %9d\n", imin, imax, log(imax-imin)/log(2), isec); }
            if (Ze < Zplane[isec])
              { imax = isec; }
            else
              { imin = isec; }
            if (ns_P != NULL) { (*ns_P)++; }
            if (imax - imin < 2) { break; }
          }
        if (use_bin || (imax - imin < 5))
          { /* Do a binary search step: */
            tosl_plane_id_t ictr = (imin + imax)/2;
            assert((imin < ictr) && (ictr < imax));
            assert(Ze != Zplane[ictr]);
            if (debug) { fprintf(stderr, "  %9d .. %9d  %5.2f bin: %9d\n", imin, imax, log(imax-imin)/log(2), ictr); }
            /* Choose sub-range and recurse: */
            if (Ze < Zplane[ictr])
              { imax = ictr; }
            else
              { imin = ictr; }
            if (ns_P != NULL) { (*ns_P)++; }
          }
      }
    assert(imax - imin == 1);
    return imax;
  }
