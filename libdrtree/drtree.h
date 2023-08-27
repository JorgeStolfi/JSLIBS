/* Drawings of asexual descent trees. */
#ifndef drtree_H
#define drtree_H
/* Last edited on 2023-06-25 18:59:05 by stolfi */

#define drtree_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

/* PURPOSE
  
  The general goal of the procedures in this library is to draw a diagram
  showing the evolution of a polulation with asexual descent.
  
  Time is assumed to be divided into discrete intervals (years, days,
  generations,...) identified by integers. Each individual in the
  population's history is also identified by a distinct index starting
  from 0. */

typedef struct drtree_node_t 
  { int32_t tbr; /* Time of birth. */
    int32_t tdt; /* Time of death. */
    int32_t par; /* Parent index, or -1. */
  } drtree_node_t; 
  /* Basic data about an individual.

    Each individual {q} in the history has a /time of birth/ {q.tbr} and a
    /time of death/ {q.tdt}, such that q.{tbr <= q.tdt}. Both may be
    negative. The /life span/ of an individual {q} is the interval of times
    {q.tbr..q.tdt}.

    Each individual {q} may have a parent {p}, whose index is {q.par}. The
    birth time of {q} must be in the range {p.tbr+1..p.tdt}}. A node {q}
    with no parent is said to be a /root/ and the notation {q.par} then
    means {NULL} by definition.

    These conditions imply that the set of all nodes in the history are
    topologically a directed forest graph, whose component are directed
    trees with a single root and one or more /leaf/ nodes (which are not
    the parent of any other node). 
  
    The /descendants/ of a node {q} comprise {q} itself
    and the descendants of all its children. 
    
    As a special case, if {q.tbr > q.tdt}, the record is a /null
    individual/ with empty life span, that should be ignored. The field
    {q.par} is then inrrelevant and should be {-1}. A null node should
    not be parent of any other node, and has no descendants (not even
    itself), by definition. */
 
#define drtree_indivs_MAX (1 << 24)
  /* Max number of individuals in diagram.  Cannot be {INT32_MAX} because
    it is used as array size. */

bool_t drtree_node_is_null(drtree_node_t *q);
  /* True iff the record {*q} describes a null individual. 
    If so, demands {q.par} to be {-1}. */

int32_t *drtree_count_children(int32_t ni, drtree_node_t dt[]);
  /* Returns a vector {nch[0..ni-1]} such that {nch[iq]} is the 
    number of individuals in {dt[0..ni-1]} that have {dt[iq]} as 
    parent.  This count will be zero for null individuals */

drtree_node_t *drtree_clip_time_range
  ( int32_t ni,
    drtree_node_t dt[],
    int32_t tMin,
    int32_t tMax
  );
  /* Given an array {dt[0..ni-1]} of basic information for individuals
    {0..ni-1}, and a time range {tMin..tMax}, extracts an array
    {dc[0..ni-1]} that describes the part of the evolution history
    comprised in the time range {tMin..tMax}, preserving the indices of
    individuals.
    
    Specifically, if the life span {q.tbr..q.tdt} of an individual {q =
    dt[iq]} is disjoint from {tMin..tMax}, then {dc[iq]} will be a null
    individual. Otherwise, the life span of {q} will be clipped to the
    range {tMin..tMax}. Moreover, if {q.tbr} is less than {tMin}, then
    {q.par} is set to {-1}, effectively turning it into a root (if it
    was not one already) -- even if part of {q}'s parent survives in {dc}. */

void drtree_check_nodes
  ( int32_t ni,
    drtree_node_t dt[],
    int32_t tMin,
    int32_t tMax
  );
  /* Checks consistency of {dt[0..ni-1]}, in particular that 
    all times are in the range {tMin..tMax}. */

#endif

