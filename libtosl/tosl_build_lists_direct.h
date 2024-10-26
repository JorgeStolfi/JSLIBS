/* Build the arc lists for topological mesh slicing, using binary+secant search. */
/* Last edited on 2024-10-07 07:08:04 by stolfi */

#ifndef tosl_build_lists_direct_H
#define tosl_build_lists_direct_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>
#include <tosl_mesh.h>

tosl_arc_id_t *tosl_build_lists_direct
  ( int32_t NP,
    tosl_coord_t Zplane[],
    tosl_mesh_t *mesh,
    int32_t debug
  );
  /* Given: the {mesh} and the list of plane Z-coordinates {Zplane[0..NP-1]}.

    Returns a array {L[0..NP-1]} of arc indices such that {L[ip]} is the
    head of the circular doubly-linked list (defined by the
    {mesh.Arc[].pred} and {mesh.Arc[].succ} links) of upward arcs in
    {mesh.Arc[0..NA-1]} that start between planes {ip-1} and {ip} and
    cross at least plane {ip} (or {-1} if there is no such arc); where
    {NA=2*NE} and {NE = mesh.NE}.
    
    For each arc {ia} that is not included in any of those lists, the
    procedure sets {Arc[ia].pred} and {Arc[ia].succ} to {ia}, makiing it
    into a singleton list.
    
    Assumes that the array Zplane[] is sorted in strictly increasing
    order. Assumes that every vertex {Z}-coordinate must be different
    from {Zplane[ip]} for all {ip}.
    
    Computes an initial plane index {ip0} by simple affine interpolation
    of each edge's {Z} between {Zplane[0]} and {Zplane[NP-1]}. Then
    adjust {ip0} by sequential search to find the correct plane index
    {ip}. The running time will be about {A*NE*NS} where {A} is some
    constant and {NS} is the max correction {|ip-ip0|} needed. Therefore,
    this method should be the best one when {Zplane[0..NP-1]} deviates
    little from the straingt line interpolation, but will lose to binary
    search when that deviation is substantial. */
    
#endif
