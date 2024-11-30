/* See {drt_lineage.h} */
/* Last edited on 2023-06-24 10:56:34 by stolfi */

#define drtree_lineage_C_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <vec.h> 
#include <affirm.h> 

#include <drtree.h>

#include <drtree_lineage.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

int32_t *drtree_lineage_collect_surviving
  ( int32_t ni,
    drtree_node_t dt[], 
    int32_t t0,
    int32_t t1
  )
  { bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }
   
    demand((ni >= 0) && (ni < drtree_indivs_MAX), "invalid {ni}");

    int32_t *fnd = (int32_t*)notnull(malloc(ni*sizeof(int32_t)), "no mem");
    
    /* Set {fnd[iq]} to -2 iff {iq} died on or after {t0} and */
    /* has surviving descendants through chidlren born after {t0}, else {-1}: */
    for (uint32_t iq = 0;  iq < ni; iq++) { fnd[iq] = -1; }
    if (debug) { fprintf(stderr, "    finding surviving lineage members\n"); }
    for (int32_t iq = ni-1; iq >= 0; iq--)
      { if ((dt[iq].tbr <= t1) && (t1 <= dt[iq].tdt))
          { if (debug) { fprintf(stderr, "      %d is a survivor\n", iq); }
            fnd[iq] = -2;
          }
        if ((fnd[iq] == -2) && (dt[iq].tbr > t0))
          { /* Has a surviving descendant but is not its own lineage. */
            if (debug) { fprintf(stderr, "      %d did not die too early and has surviving descs\n", iq); }
            int32_t ip = dt[iq].par;
            if (ip != -1)
              { demand((ip >= 0) && (ip < iq), "bad parent index");
                if (dt[ip].tdt >= t0)
                  { if (debug) { fprintf(stderr, "      propagating to parent %d\n", ip); }
                    fnd[ip] = -2;
                  }
              }
          }
      }
      
    /* Now scan forward assigning founder indices: */
    if (debug) { fprintf(stderr, "    finding lineage founders\n"); }
    for (uint32_t iq = 0;  iq < ni; iq++) 
      { if (fnd[iq] == -2)
          { if (debug) { fprintf(stderr, "      %d did not die too early and has surviving descs\n", iq); }
            assert(dt[iq].tdt >= t0);
            fnd[iq] = -1; /* By default, unless... */
            if (dt[iq].tbr <= t0)
              { /* {iq} is a founder. */
                if (debug) { fprintf(stderr, "        %d is a founder\n", iq); }
                fnd[iq] = iq;
              }
            else
              { /* {iq} may be descendant of a founder: */
                int32_t ip = dt[iq].par;
                if (ip != -1)
                  { assert((ip >= 0) && (ip < iq));
                    if (fnd[ip] >= 0)
                      { if (debug) { fprintf(stderr, "        %d, like its parent %d, belongs to %d's lineage\n", iq, ip, fnd[ip]); }
                        fnd[iq] = fnd[ip];
                      }
                    else
                      { if (debug) { fprintf(stderr, "        %d's parent %d has no lineage\n", iq, ip); }
                      }
                  }
                else
                  { if (debug) { fprintf(stderr, "        %d has no parent and is not founder\n", iq); } }
              }
          }
        else
          { assert(fnd[iq] == -1); }
      }
      
    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
    return fnd;
  }
