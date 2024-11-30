/* See enum_orbits.h. */
/* Last edited on 2024-11-23 06:05:38 by stolfi */

#define enum_orbits_C_copyright \
  "Copyright © 1996, 2006 Institute of Computing, Unicamp."

/* AUTHORS

  These enumeration procedures were originally developed as the
  Modula-3 file {Oct.m3} by J. Stolfi and Rober M. Rosi in 1993. The
  latter was converted to C by J. Stolfi in 1996 and was substantially
  revised by him in January 2007. */

#include <stdint.h>
#include <malloc.h>
#include <assert.h>

#include <jswsize.h>
#include <ulist.h>
#include <ref.h>
#include <bool.h>
#include <vec.h>

#include <enum_orbits.h>

bool_t enum_cycle(ref_t p, enum_step_t step, enum_visit_t visit, ref_vec_t *vP)
  { uint32_t nP = 0; /* Visited places are {vP->e[0..nP-1]}, if {vP != NULL}. */ 
    ref_t q = p;
    bool_t stop = FALSE;
    do 
      { stop = ( visit == NULL ? FALSE : visit(q));
        if (vP != NULL) { ref_vec_expand(vP, nP); vP->e[nP] = q; nP++; }
        q = step(q);
      } 
    while ((q != p) && (! stop));
    if (vP != NULL) { ref_vec_trim(vP, nP); }
    return stop;
  }

#define INIT_QUEUE_SIZE 1024

bool_t enum_items
  ( ref_vec_t root, 
    uint32_t ns, 
    enum_step_t *step[], 
    enum_visit_t *visit[], 
    ref_vec_t *vP
  )
  {
    ref_vec_t Q = (vP != NULL ? *vP : ref_vec_new(INIT_QUEUE_SIZE));
    uint32_t nQ = (vP != NULL ? vP->ne : 0);
      /* The places visited so far (including pre-existing entries
        in {vP}, f any) are {Q[0..nQ-1]}, in stratified order. */

    auto void init_marking(void);
      /* Initializes the tables used by {Marked} and {Mark}. */
    
    auto bool_t check_and_set_mark(ref_t p);
      /* If {p} is not marked as visited, marks {p} and returns FALSE.
        If {p} is already marked, returns TRUE. */
    
    auto void done_marking(void);
      /* Clears any permanent marks and reclaims storage used by
        {is_marked} and {mark}. */
    
    bool_t stop = FALSE; /* Set to TRUE if any {visit[k]} returns TRUE. */

    auto void enum_one(ref_t p, uint32_t n);
      /* If {p} is already marked, do nothing. If {p} is unmarked, the
        procedure assumes that {p} was reached through {step[n]} (or,
        if {n==nst}, that {p} is a root). The procedure then marks and
        appends to {Q} the node {p} and every node {q} that can be
        reached from {p} by chains of one or more steps in
        {step[0..n-1]}, visiting each {q} with {visit[k]}, where {k}
        is the index of the last step used to reach {q}. Stops as soon
        as a {visit[k]} procedure returns TRUE. */
    
    auto void close_orbit(uint32_t n);
      /* Finds all unmarked items that can be reached from item
        {Q[nQ-1]} by chains of steps in {step[0..n-1]}; marks every
        such item {q}, appends {q} to {Q[0..nQ-1]}, and visits {q}
        with {visit[k]} where {k} is the index of the last step used
        to reach {q}. Stops as soon as a {visit[k]} procedure returns
        TRUE. */
    
    /* Initialize the mark tables: */
    init_marking();
    
    /* Enumerate all roots: */
    for (uint32_t i = 0;  (i < root.ne) && (! stop); i++)
      { enum_one(root.e[i], ns); }

    /* Release the marking table storage: */
    done_marking();
    
    /* Return {Q}, or reclaim its space: */
    if (vP == NULL)
      { free(Q.e); }
    else
      { ref_vec_trim(&Q, nQ); (*vP) = Q; }
      
    /* Return the aborted-enum indicator: */
    return stop;

    /* MAIN PROCEDURE */
    
    void enum_one(ref_t p, uint32_t n)
      { if (! check_and_set_mark(p))
          { /* New item, stack it: */
            ref_vec_expand(&Q, nQ); 
            Q.e[nQ] = p; nQ++;
            /* Visit it: */
            if (visit[n] != NULL) 
              { stop = visit[n](p); if (stop) { return; } } 
            /* Gather its orbit under {step[0..n-1]}: */
            close_orbit(n);
          }
      }
    
    void close_orbit(uint32_t n)
      { if (n > 0)
          { assert(nQ > 0);
            uint32_t nX = nQ - 1;
            close_orbit(n-1);
            while((nX < nQ) && (! stop))
              { ref_t q = Q.e[nX]; nX++;
                q = step[n-1](q);
                enum_one(q, n-1);
              }
          }
      }
             
    /* MARKING PROCEDURES */
    
    ulist_t *M;
      /* The set {M} contains the addresses visited so far. */
    
    void init_marking(void)
      { M = ulist_new(INIT_QUEUE_SIZE); }
    
    bool_t check_and_set_mark(ref_t p)
      { uaddress_t up = (uaddress_t)p;
        if (ulist_has(M, up))
          { return TRUE; }
        else
          { (void)ulist_insert_last(M, up);
            uint32_t ctM = ulist_count(M);
            uint32_t szM = ulist_capacity(M);
            if (2*ctM > szM) { ulist_resize(M, 2*szM+1, NULL); }
            return FALSE;
          }
      }
    
    void done_marking(void)
      { ulist_free(M); }
  }

bool_t enum_orbits
  ( ref_vec_t root,
    uint32_t ni, 
    enum_step_t *istep[], 
    uint32_t no, 
    enum_step_t *ostep[], 
    enum_visit_t *visit,
    ref_vec_t *vP
  )
  {
    uint32_t nvis; /* Counts arcs added to {vP}. */
    
    auto bool_t do_visit(ref_t p);
    /* A visit-func to be called when starting a new root or 
       a new {istep}-orbit.  It appends {p} to {vP} at position
       {nvis}, calls the client-given visit-func, and increments
       {nvis}. */
    
    /* Create the step-funcs and visit-funcs for enumeration: */
    uint32_t ns = ni + no;
    enum_step_t *stp[ns];
    enum_visit_t *vis[ns+1];
    for (uint32_t i = 0;  i < ns; i++)
      { stp[i] = (i < ni ? istep[i] : ostep[i - ni]);
        vis[i] = (i < ni ? NULL : &do_visit);
      }
    vis[ns] = &do_visit;
    
    nvis = 0;
    bool_t res = enum_items(root, 3, stp, vis, NULL);
    if (vP != NULL) { ref_vec_trim(vP, nvis); }
    return res;
    
    /* IMPEMENTATION OF LOCAL FUNCTIONS */
    
    bool_t do_visit(ref_t p)
      { if (vP != NULL) 
          { ref_vec_expand(vP, nvis);
            vP->e[nvis] = p;
          }
        nvis++;
        return (visit == NULL ? FALSE : visit(p));
      }
  }

