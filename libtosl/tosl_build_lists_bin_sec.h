/* Build the arc lists for topological mesh slicing, using binary+secant search. */
/* Last edited on 2024-10-06 16:47:09 by stolfi */

#ifndef tosl_build_lists_bin_sec_H
#define tosl_build_lists_bin_sec_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>
#include <tosl_mesh.h>

tosl_arc_id_t *tosl_build_lists_bin_sec
  ( int32_t NP,
    tosl_coord_t Zplane[],
    tosl_mesh_t *mesh,
    int32_t use_bin,
    int32_t use_sec,
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
    order. The array {Zedge} may be in arbitrary order. Assumes that
    {Zedge[ie]} is different from {Zplane[ip]} for all {ie,ip}.
    
    If the booleans {use_bin} and {use_sec} are both true (normal case)
    uses a combination of binary and secant search. The run time should
    be almost {\O(NE+NP)} if the spacing of the planes has limited
    variation.
    
    If only {use_bin} is true, uses pure binary search, which should
    take time {\O(NP + NE\log NP)} for any plane spacing. If only
    {use_sec} is true, uses pure secant search, which may take
    {\THETA(NP*NE)} if the plane spacing is perverse enough. */
    
#endif
