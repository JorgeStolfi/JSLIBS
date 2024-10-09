/* Manipulating the linked lists of arcs in a {tosl_mesh_t}. */
/* Last edited on 2024-10-08 22:55:09 by stolfi  */

#ifndef tosl_arc_list_H
#define tosl_arc_list_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <tosl.h>
#include <tosl_mesh.h>
  
/* LINKED LIST MANIPULATION

  The following procedures handle circular doubly-linked lists of 
  {tosl_arc_t} records in the array {mesh.Arc[]}, formed by the
  links {.pred} and {.succ} links in {mesh.Arc[]}.*/

int32_t tosl_arc_list_len(tosl_arc_id_t L, tosl_mesh_t *mesh);
  /* Returns the number of elements in the lust that starts at {Arc[L]}. */

tosl_arc_id_t tosl_arc_list_pop(tosl_arc_id_t *L_P, tosl_mesh_t *mesh);
  /* Removes an element from a list.
  
    Expects {ia=*L_P} to be the index into {mesh.Arc[]} of one element of a
    non-empty arc list. Removes element {mesh.Arc[ia]} from that list and sets
    {*L_P} to the index of an element of what remains of the list
    (possibly {-1}, if that was the last element). Then makes {mesh.Arc[ia]}
    into a singleton list and returns the index {ia}. */
    
void tosl_arc_list_remove(tosl_arc_id_t *L_P, tosl_arc_id_t ia, tosl_mesh_t *mesh);
  /* Removes arc {ia} from its list.
  
    Expects {*L_P} to be the index into {mesh.Arc[]} of one element of a
    non-empty arc list, and {ia} to be the index of some element of that 
    list. Removes element {mesh.Arc[ia]} from that list and makes it into a 
    singleton list.  If {ia} was {*L_P},
    sets {*L_P} to the index of some other element of that list
    (possibly {-1}, if that was the last element). */ 
    
void tosl_arc_list_add(tosl_arc_id_t *L_P, tosl_arc_id_t ia, tosl_mesh_t *mesh); 
  /* Adds an element to a list.
  
    Expects {mesh.Arc[ia]} to be the single element of a list, and {*L_P}
    to be the index of an element of another list {L}, or {-1} to mean
    that {L} is the empty list. The list {L} must not contain {Arc[ia]}.
    If {L} is not empty, inserts {mesh.Arc[ia]} into that list and leaves
    {*L_P} unchanged. If {L} is empty, just sets {*L_P} to {ia}. */

tosl_arc_id_t tosl_arc_list_merge
  ( tosl_arc_id_t *L0_P,
    tosl_arc_id_t *L1_P,
    tosl_coord_t Zp,
    tosl_mesh_t *mesh
  );
  /* Given the head indices {*L0_P} and {*L1_P} of two circular
    doubly-linked lists of arcs {L0} and {L1}, combines them into a
    single list {L}, and returns the index of its head element.
    
    The two lists are assumed to be distinct, hence disjoint. Every arc
    in either list is expected to be directed upwards and to start below
    the plane with {Z}-coordinate {Zp}. Elements of {L0} that end below
    {Zp} are left out. Elements of {L1} are assumed to end above {Zp}.
    Thus the final list {L} will have only upward arcs that cross that
    plane.
    
    Either input index {*L0_P} or {*L1_P} may be {-1} to signify that it
    is empty. If both are empty, returns {-1}. In any case, o return
    {*L0_P} and {*L1_P} will be set to {-1}. */

#endif
