/* Build the arc lists for topological mesh slicing, using binary+secant search. */
/* Last edited on 2024-10-07 07:07:33 by stolfi */

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
    order.  Every vertex {Z}-coordinate must be
    different from {Zplane[ip]} for all {ip}.
    
    If the booleans {use_bin} and {use_sec} are both true (normal case)
    uses a combination of binary and secant search. The run time should
    be almost {A1*NE}, for some constant {A1}, if the spacing of the
    planes has limited variation, and {A2*NE*log(NP)} in the worst case.
    
    If only {use_bin} is true, uses pure binary search, which, for any
    plane spacing, should take time {A3*NE*log(NP)}, for some constants
    {A3,B} with {A3<A2}. If only {use_sec} is true, uses pure secant
    search, which may take {A4*NP*NE)}, for some constant {A3}, if the
    plane spacing is perverse enough. */
    
#endif
