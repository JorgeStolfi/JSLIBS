/* See {drtree_test.h}.  */
/* Last edited on 2023-06-24 10:59:16 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <cmp.h>
#include <affirm.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <argparser.h>
#include <vec.h>  

#include <drtree.h>

#include <drtree_test.h>
  
void drtree_test_create_individuals
  ( int32_t ni,       /* Number of individuals in population. */
    int32_t nRoots,   /* Number of root individuals. */
    int32_t tStart,   /* Nominal start time of evolution. */
    int32_t tStop,    /* Nominal final time of evolution. */
    bool_t orphans,   /* If false, no roots allowed after {tStart}. */
    int32_t ageMax,   /* Max individual age. */
    int32_t nchMax,   /* Max num of children. */
    drtree_node_t dt[]
  )
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }
    
    demand((ni >= 0) && (ni < drtree_indivs_MAX), "invalid {ni}");

    int32_t nch[ni]; /* Child counts. */
    int32_t nr = 0; /* Number of roots actually created. */
    demand (ni == ni, "inconsistent {ni}");
    int32_t ni_act = 0; /* number of indivs actually created. */
    
    /* Range of all times seen: */
    int32_t tMin = INT32_MAX;
    int32_t tMax = INT32_MIN;

    auto int32_t choose_parent(int32_t iq);
      /* Chooses a valid parent for node {iq}, possibly {-1}. If {-1},
        then there are not too many roots yet, and {iq} should be a root.
        Othwerwise the result will have a long enough life span and not too
        many children already. */

    auto bool_t is_good_parent(int32_t jq, int32_t iq);
      /* True if {jq} (which must be in {0..iq-1}) is a viable parent for {iq}, namely 
        has a long enough life span and not too many children already.
        Does not check limits on number of roots. */

    while (ni_act < ni)
      { /* Index of next node: */
        int32_t iq = ni_act;
        drtree_node_t *q = &(dt[iq]);
        /* Choose the age at death {L} */
        int32_t L = int32_abrandom(0, ageMax);
        /* Decide on the parent {pari}: */
        int32_t ip =  choose_parent(iq);
        
        /* Choose the birth time: */
        int32_t tbri;     
        if (ip == -1)
          { assert(nr < nRoots);
            int32_t tbrMin = tStart - L; /* Earliest birth time allowed. */
            int32_t tbrMax = (orphans ? tStop+1 : tStart); /* Earliest birth time allowed. */
            tbri = int32_abrandom(tbrMin, tbrMax);
          }
        else
          { assert((ip >= 0) && (ip < iq));
            drtree_node_t *p = &(dt[ip]);
            assert(p->tdt > p->tbr);
            assert(p->tdt - p->tbr >= nchMax);
            assert(nch[ip] < nchMax);
            tbri = int32_abrandom (p->tbr+1, p->tdt);
          }
        /* Choose the death time: */
        int32_t tdti = tbri + L;
        /* Save: */
        q->tbr = tbri; if (tbri < tMin) { tMin = tbri; }
        q->tdt = tdti; if (tdti > tMax) { tMax = tdti; }
        q->par = ip;
        nch[iq] = 0; /* to be incremented. */
        fprintf(stderr, "    node %d life span %d {%d..%d} parent %d", iq, L, tbri, tdti, ip);
        if (ip == -1)
          { nr++; }
        else
          { nch[ip]++;
            fprintf(stderr, " nch = %d", nch[ip]);
          }
        fprintf(stderr, "\n");
        ni_act++;
      }

    if (debug) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
    return;
    
    int32_t choose_parent(int32_t iq)
      { 
        for (uint32_t it = 0;  it < 1000; it++)
          { /* Choose node {s} earlier than {iq}, favoring close to {iq}: */
            int32_t sMin = (nr < nRoots ? nr - nRoots : 0);
            int32_t sMax = iq - 1;
            int32_t sMid1 = int32_abrandom(sMin, sMax);
            int32_t sMid2 = int32_abrandom(sMid1, sMax);
            int32_t sMid3 = int32_abrandom(sMid2, sMax);
            int32_t sMid4 = int32_abrandom(sMid3, sMax);
            int32_t s = int32_abrandom(sMid4, sMax);
            if (debug) { fprintf(stderr, "      s = %d", s); }
            /* Choose a viable candidate parent {ip}, either {s}'s parent, {s}, or {-1}: */
            int32_t ip;
            if (s < 0)
              { ip = -1; }
            else 
              { /* Try using {s} as sibling, if not as parent, if not retry: */
                ip = dt[s].par;
                if (debug) { fprintf(stderr, " ip = %d", ip); }
                if ((ip == -1) || (! is_good_parent(ip, iq))) 
                  { ip = s; 
                    if (debug) { fprintf(stderr, " ip = %d", ip); }
                    if (! is_good_parent(ip, iq)) 
                      { ip = -1;
                        if (debug) { fprintf(stderr, " ip = %d", ip); }
                      }
                  }
              }
            if (debug) { fprintf(stderr, "\n"); }
            /* If within the limits, that is it: */
            if ((ip != -1) || (nr < nRoots)) { return ip; }
          }
        demand(FALSE, "failed to choose a parent");
      }
      
    bool_t is_good_parent(int32_t jq, int32_t iq)
      { assert((jq >= 0) && (jq < iq));
        int32_t Lj = dt[jq].tdt - dt[jq].tbr + 1; /* Number of times in life span of {jq}. */
        if (Lj <= nchMax)
          { /* Life span of {jq} is not long enough: */
            if (debug) { fprintf(stderr, " life span (%d) = %d times", jq, Lj); }
            return FALSE;
          }
        else if (nch[jq] >= nchMax) 
          { /* Node {jq} already has too many children: */
            if (debug) { fprintf(stderr, " nch (%d) = %d", jq, nch[jq]); }
            return FALSE;
          }
        else
          { return TRUE; }
      }
    
    drtree_check_nodes(ni, dt, tMin, tMax);

  }

drtree_test_options_t *drtree_test_parse_options(argparser_t *pp)
  { 
    drtree_test_options_t *o = notnull(malloc(sizeof(drtree_test_options_t)), "no mem");

    argparser_get_keyword(pp, "-nIndivs");
    o->nIndivs = (int32_t)argparser_get_next_int(pp, 1, drtree_test_num_indivs_MAX);

    argparser_get_keyword(pp, "-nRoots");
    o->nRoots = (int32_t)argparser_get_next_int(pp, 1, o->nIndivs);

    argparser_get_keyword(pp, "-tStart");
    o->tStart = (int32_t)argparser_get_next_int(pp, -100000, +3000);

    argparser_get_keyword(pp, "-tStop");
    int32_t tStopMax = o->tStart + drtree_test_num_times_MAX - 1;
    o->tStop = (int32_t)argparser_get_next_int(pp, o->tStart+2, tStopMax);

    argparser_get_keyword(pp, "-orphans");
    o->orphans = argparser_get_next_bool(pp);

    argparser_get_keyword(pp, "-tRef");
    o->tRef = (int32_t)argparser_get_next_int(pp, o->tStart,o->tStop);

    argparser_get_keyword(pp, "-ageMax");
    o->ageMax = (int32_t)argparser_get_next_int(pp, 1, drtree_test_num_times_MAX);

    argparser_get_keyword(pp, "-nchMax");
    o->nchMax = (int32_t)argparser_get_next_int(pp, 1, drtree_test_num_children_MAX);

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);
    return o;
  }
