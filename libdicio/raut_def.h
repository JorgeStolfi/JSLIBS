#ifndef raut_def_H
#define raut_def_H

#define raut_def_DESC "Internal representation of {raut_t}"

/* Last edited on 2009-10-28 23:42:37 by stolfi */
/*©*/

#include <stdint.h>

#include <bool.h>
#include <rdag.h>
#include <raut.h>

struct raut_t
  { 
    char *doc;              /* A multi-line comment (NULL, empty, or ends with '\n\000'). */
    raut_state_t root; /* The current root state. */
    rdag_t *D;              /* The underlying dag. */
    
    /* Auxiliary data for finding prefixes: */
  };
  
/*
  UNSAFE PROCEDURES: The following procedures may break the invariants
  of the automaton, and should be used only by people in the know.
  
  REPRESENTATION: An automaton {A} is represented here by a dag {D}
  with some restrictions and some extra info and intepretation.
  
  The restrictions are:
  
  * the output symbols of {D} are the bits 0 and 1
  
  * Any proper node {s} of {D} whose p-link is null has its o-mark equal to one. 
  
  The extra info is 
  
  * A distinguished /root state/ {u_root = (o_root,t_root)}.
  
  The extra interpretations are:
  
  * The input symbols of the automaton {A} are those of {D}.
  
  * A state of the automaton {A} is a pair {u == (o,s)} consisting of an
    output symbol {o} of {D} and a node {s} of {D}.
  
  * The arcs out of a state {u == (o,s)} f {A} are all the 
    triplets  {(i,u') == (i,(o',s'))} such that there exists a proper 
    subnode {t} of {s} with i-mark {i}, o-mark {o'}, and p-link {s'}. 
    
  Note that the the outgoing arcs of a state {u=(o,s)} depend only on
  the node {s} and not on the but {o}. Note also that the f-links
  traverse the arc list in reverse order. More precisely, if state
  {u=(o,s)} is neither void nor unit, its *last* outgoing arc consists
  of the i-mark, o-mark, and p-link of the node {s} itself. Moreover,
  if {t} is a subnode of {s}, the arcs out of the state {(0,t)} or
  {(1,t)} are an initial prefix of those out of {u}.

  Any path {P = (u[0],a[1],u[1],...,a[n],u[n])} in the automaton {A}
  corresponds to a path {P' = (o[0],t[0],s[1],i[1],o[1],t[1],...,s[n],i[n],o[n],t[n])} 
  in the dag {D}.  Namely if {P} starts at the state {u[0]=(o',s')}, then {P'} starts
  with {o[0]=o'} and {t[0]=t'}.  If the {k}th arc {a[k]} of {P} is {(i[k],(o[k],t[k]))}, then
  the {k}th step of {P'} is {(s[k],i[k],o[k],t[k])}, where {s[k]} is 
  the unique subnode of {t[k-1]} that has i-mark {i[k]}.  */
  
rdag_t *raut_dag_get(raut_t *A);
  /* The dag underlying the automaton {A}.  Note that 
    modifing {D} directly may break the invariants of {A}. */

#endif

//        (*************************************************************************)
//        (* PREFIX/SUFFIX DATA                                                    *)
//        (*************************************************************************)
//  
//        (*
//          The folowing fields are auxiliary tables used by methods that deal
//          with prefixes and suffixes.
//  
//          The "nSuffs" "nSufLetters" arrays are always allocated, big enough, 
//          and up-to-date.
//  
//          The prefix-related fields ("rdag", "rev", "dir", "nPrefs", and 
//          "nPrefLetters") are (re)computed only when needed (and if "aut.root" 
//          is not raut_state_NULL). If "prefRoot" is raut_state_NULL, then those fields are
//          garbage; if "prefRoot" is a proper state, those fields are valid
//          provided that the current root is "prefRoot".
//  
//          The reverse DAG "rdag" has one arc from rev[t] to rev[s]
//          whenever the automaton has a proper arc from "s" to "t".
//          The virtual arcs of "aut" that go to raut_state_NULL, including
//          the invisible arcs labelled with NullLetter that denote final states
//          of "aut", have no correspondents in "rdag". On the other hand, "rdag"
//          has one arc labelled with "NullLetter" going from rev[prefRoot]
//          to raut_state_NULL. *)
//  
//        nSuffs: REF ARRAY OF uint32_t;       (* Number of suffix strings per state *)
//        nSuffLetters: REF ARRAY OF uint32_t; (* Total length of suffix strings per state *)
//  
//        prefRoot: raut_state_t;     (* raut_root used to compute "rdag, rev, dir, nPrefs" *)
//  
//        rdag: rdag_t;                    (* The DAG with reverse transitions *)
//        rev: REF ARRAY OF DAG.raut_state_t;    (* Maps states of self to states of "rdag" *)
//        dir: REF ARRAY OF raut_state_t;        (* Inverse of "rev", or raut_state_NULL if undef *)
//        nPrefs: REF ARRAY OF uint32_t;       (* Number of prefix strings for each state *)
//        nPrefLetters: REF ARRAY OF uint32_t; (* Total length of prefix strings per state *)

