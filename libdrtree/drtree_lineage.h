/* Identifying lineages in asexual descent trees. */
#ifndef drtree_lineage_H
#define drtree_lineage_H
/* Last edited on 2023-06-15 13:06:50 by stolfi */

#define drtree_lineage_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <vec.h> 

#include <drtree.h> 

/* 
  The following concepts appply to a tree :
  
  LINEAGE
  
  The /lineage/ of an individual {q} is the set {L(q)} of all its
  descendants, including {q} itself.
   
  Note that two distinct individuals born on the same year will have
  disjoing lineages.
 
  LINEAGE AS OF
  
  If an individual {q} is alive at some time {t0} (that is, {q.tbr <=
  t0 <= q.tdt}), its /lineage as of/ {t0} is the set{L(q|t0)}
  consisting of itself plus {L(r)} for all children {r} born strictly
  after {t0}.  
  
  If {q} is not alive at time {t0}, then {L(q|t0)} is empty by convention.
  
  Note that {L(q|t0)} is a subset of {L(q)}, but is disjoing from the lineage 
  {L(c)} of any child {c} of {q} that is born on or before time {t0}.
  
  The individual {q} is said to be the /founder/ of the lineage {L(q|t0)}. 
  
  SURVIVING LINEAGES
  
  We say that a lineage {L(q)} is /surviving at/ some time {t1} if it has some
  individual {r} who is alive at {t1} - that is, with {r.tbr <= t1 <= r.tdt}.
  
  For times {t0,t1} with {t0 < t1}, we define {L(q|t0|t1)}, the
  /lineage/ of {q} /as of/ {t0} /surviving to/ {t1}, as the set of those
  individuals {r} in {L(q|t0)} such that {L(r)} is surviving at {t1}.
  If {L(q,t0)} is empty, or there is no such {r}, then {L(q|t0|t1)} 
  is empty.
  
  FOUNDERS
  
  The individual {q} is said to be the /founder/ of the lineage {L(q)},
  {L(q|t0)}, or {L(q|t0|t1)} if the set is not empty.  */

int32_t *drtree_lineage_collect_surviving
  ( int32_t ni,
    drtree_node_t dt[], 
    int32_t t0,
    int32_t t1
  );
  /* Identifies each individual {q=dt[iq]} in {dt[0..ni-1]} which is
    the founder of a non-empty lineage {L(q|t0|t1)}. 
    
    Returns a list {fnd[0..ni-1]} such that {fnd[ir]} is {iq} iff
    individual {r = dt[ir]} belongs to the lineage {L(q|t0|t1)}
    where {q=dt[iq]}. Otherwise {fnd[ir]} is set to {-1}.
    
    In particular, {fnd[ir]} is {-1} if individual {r=dt[ir]} died before
    {t0}, has no ancestor born on or before {t0}, or has no descendant
    who is alive at time {t1}. */

#endif

