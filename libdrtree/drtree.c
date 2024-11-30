/* See {drtree.h} */
/* Last edited on 2023-06-24 11:04:47 by stolfi */

#define drtree_C_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h> 
#include <vec.h> 
#include <jsmath.h> 

#include <drtree.h> 

/* INTERNAL PROTOTYPES */
    
/* IMPLEMENTATIONS */

bool_t drtree_node_is_null(drtree_node_t *q)
  { bool_t res = (q->tbr > q->tdt);
    if (res) { demand(q->par == -1, "invalid {par} in null node"); }
    return res;
  }

int32_t *drtree_count_children(int32_t ni, drtree_node_t dt[])
  { 
    demand((ni >= 0) && (ni < drtree_indivs_MAX), "invalid {ni}");
    
    int32_t *nch = (int32_t *)notnull(malloc(ni*sizeof(int32_t)), "no mem");
    for (uint32_t iq = 0;  iq < ni; iq++) { nch[iq] = 0; }
    
    /* Scan in reverse chrono order, defining the {nch} fields: */
    for (int32_t iq = ni-1; iq >= 0; iq--)
      { drtree_node_t *q = &(dt[iq]);
        if (! drtree_node_is_null(q))
          { int32_t ip = q->par;
            if (ip != -1)
              { /* Node {iq} is a child of {ip}. */
                /* Consistency checks: */
                assert((ip >= 0) && (ip < iq));
                drtree_node_t *p = &(dt[ip]);
                assert(! drtree_node_is_null(p)); /* Null indiv can't be parent. */
                assert(p->tbr <= q->tbr); /* Can't be parent before birth. */
                assert(p->tdt >= q->tbr); /* Can't be parent after death. */
                nch[ip]++;
              }
          }
      }

    return nch;
  }

drtree_node_t *drtree_clip_time_range
  ( int32_t ni,
    drtree_node_t dt[],
    int32_t tMin,
    int32_t tMax
  )
  {
    demand((ni >= 0) && (ni < drtree_indivs_MAX), "invalid {ni}");
    demand(tMin < tMax, "invalid {tMin..tMax}");

    drtree_node_t *dc = (drtree_node_t*)notnull(malloc(ni*sizeof(drtree_node_t)), "no mem");
    for (uint32_t iq = 0;  iq < ni; iq++)
      { drtree_node_t *q_in = &(dt[iq]);
        drtree_node_t *q_ot = &(dc[iq]);
        if (drtree_node_is_null(q_in))
          { (*q_ot) = (*q_in); }
        else
          { q_ot->tdt = (q_in->tdt > tMax ? tMax : q_in->tdt);
            q_ot->tbr = (q_in->tbr < tMin ? tMin : q_in->tbr);
            if (q_ot->tbr > q_ot->tdt)
              { /* {q_in} is entirely outside {tMin..tMax}, becomes null node: */
                q_ot->par = -1;
              }
            else if (q_in->tbr < tMin)
              { /* {q_in} birth is before {tMin}, becomes root: */
                q_ot->par = -1;
              }
            else
              { /* Birth of {q} is in range, preserve parent: */
                int32_t ip = q_in->par;
                if (ip != -1)
                  { /* Paranoia, check validity of parent: */
                    assert((ip >= 0) && (ip < iq));
                    drtree_node_t *p_in = &(dt[ip]);
                    assert(! drtree_node_is_null(p_in));
                    assert((p_in->tbr <= q_in->tbr) && (q_in->tbr <= p_in->tdt));
                  }
                q_ot->par = ip;
              }
          }
      }
    return dc;
  }
            
void drtree_check_nodes
  ( int32_t ni,
    drtree_node_t dt[],
    int32_t tMin,
    int32_t tMax
  )
  {
    demand((ni >= 0) && (ni < drtree_indivs_MAX), "invalid {ni}");
    demand(tMin < tMax, "invalid {tMin..tMax}");
    for (uint32_t iq = 0;  iq < ni; iq++)
      { drtree_node_t *q = &(dt[iq]);
        if (drtree_node_is_null(q))
          { /*Life span must be empty: */
            assert(q->tbr > q->tdt);
            /* Must not have a parent: */
            assert(q->par == -1);
          }
        else
          { /* Life span must be non-empty and in {tMin..tMax}: */
            assert((tMin <= q->tbr) && (q->tbr <= q->tdt) && (q->tdt <= tMax));
            /* Check its parent: */
            int32_t ip = q->par;
            if (ip != -1)
              { /* Parent must come before children: */
                assert((0 <= ip) && (ip < iq));
                drtree_node_t *p = &(dt[ip]);
                /* Parent can't be a null node: */
                assert(! drtree_node_is_null(p));
                /* Child must be born during parent's life span: */
                assert((p->tbr <= q->tbr) && (q->tbr <= p->tdt));
              }
          }
      }
  }
          
