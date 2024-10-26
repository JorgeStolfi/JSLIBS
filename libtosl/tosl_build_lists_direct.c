/* See {tosl_build_lists_direct.h}. */
/* Last edited on 2024-10-07 14:58:00 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <tosl.h>
#include <tosl_mesh.h>
#include <tosl_arc_list.h>
#include <tosl_build_lists_direct.h>

tosl_plane_id_t tosl_build_lists_direct_search
  ( tosl_coord_t Ze,
    int32_t NP,
    tosl_coord_t Zplane[],
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
    
    If {ns_P} is not {NULL}, increments it with the number of sequential steps 
    performed. */

tosl_arc_id_t *tosl_build_lists_direct
  ( int32_t NP,
    tosl_coord_t Zplane[],
    tosl_mesh_t *mesh,
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
            tosl_plane_id_t ip = tosl_build_lists_direct_search(Zorg, NP, Zplane, debug, ns_P);
            /* Paranoia, checking result of plane search: */
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
      { double nspe = ((double)ns)/((double)mesh->NE);
        fprintf(stderr, "performed %.2f sequential steps per edge\n", nspe);
      }
    return L;
  }

tosl_plane_id_t tosl_build_lists_direct_search
  ( tosl_coord_t Ze,
    int32_t NP,
    tosl_coord_t Zplane[],
    int32_t debug,
    int32_t *ns_P
  )
  {
    if (debug) { fprintf(stderr, "  Ze = %+10d\n", Ze); } 
    tosl_coord_t Z0 = Zplane[0], Z1 = Zplane[NP-1];
    if (Ze < Z0) { return 0; }
    if (Ze > Z1) { return NP; }
    assert((Z0 < Ze) && (Ze < Z1)); /* No plane-vertex collision. */
        
    /* Use affine interpolation to guess {ip}: */
    double dzedge = Ze - Z0;
    double dzplane = Z1 - Z0;
    tosl_plane_id_t ip = (tosl_plane_id_t)floor((dzedge*NP)/dzplane);
    assert((0 <= ip) && (ip < NP));

    /* Adjust {ip} by sequential search: */
    while ((ip < NP) && (Ze > Zplane[ip])) { ip++; if (ns_P != NULL) { (*ns_P)++; } }
    while ((ip > 0) && (Ze < Zplane[ip-1])) { ip--; if (ns_P != NULL) { (*ns_P)++; } }
    return ip;
  }
