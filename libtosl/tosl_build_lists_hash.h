/* Build the arc lists for topological mesh slicing, using a hash table. */
/* Last edited on 2024-10-07 07:00:54 by stolfi */

#ifndef tosl_build_lists_hash_H
#define tosl_build_lists_hash_H

#define _GNU_SOURCE
#include <stdint.h>

#include <tosl.h>
#include <tosl_mesh.h>

tosl_arc_id_t *tosl_build_lists_hash
  ( int32_t NP,
    tosl_coord_t Zplane[],
    tosl_mesh_t *mesh,
    int32_t debug
  );
  /* Given: the {mesh} and the list of plane Z-coordinates {Zplane[0..NP-1]}.
    
    Returns a table {L[0..NP-1]} of arc indices such that {L[ip]} is the
    head of the circular doubly-linked list (defined  by the {Arc[].pred} and
    {Arc[].succ} links) of upward arcs that start between planes {ip-1} and
    {ip} and cross at least plane {ip} (or {-1} if there is no such
    arc).
    
    For each arc {ia} that is not included in any of those lists, the 
    procedure sets {Arc[ia].pred} and {Arc[ia].succ} to {ia},
    makiing it into a singleton list.
    
    Assumes that the array {Zplane[]} is sorted in strictly increasing
    order. Every vertex {Z}-coordinate must be different from {Zplane[ip]}
    for all {ip}.
    
    Uses a hash table method that should take time {A*NE*NC + B*NP},
    for some constants {A,B}, where {NC} is the max number of 
    {Zplane[]} values that hash to the same table entry.  This number
    should be inversely proportional to the ratio between the average spacing
    {(Zplane[NP-1]-Zplane[0])/(NP-1)} and the minimum spacing between planes.
    
    The boolean {debug} causes diagnostic printouts if true. */

#endif
