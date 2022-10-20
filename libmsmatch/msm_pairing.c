/* See msm_pairing.h */
/* Last edited on 2022-10-20 07:59:10 by stolfi */

#define msm_pairing_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <sign.h>
#include <fget.h>
#include <nget.h>
#include <vec.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <rmxn.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>
#include <msm_rung.h>
#include <msm_pairing.h>

/* INTERNAL PAIRING REPRESENTATION */

struct msm_pairing_t 
  { int32_t ng;             /* Number of rungs in pairing. */
    msm_rung_t gini;    /* Rung number 0. */
    msm_rung_vec_t gv;  /* Rungs of pairing, or empty if perfect. */
  };
  /* The number of rungs in the pairing is {ng}, which is 
    non-negative.
    
    If {ng} is zero, the pairing is empty (has no rungs).
    In that case {gini} is irrelevant, usually {(0,0)}.
    
    If {ng} is positive but {gv.ne} is 0, the pairing is implicit and
    perfect -- namely, {pr[i] = gini + (i,i)} for all {i} in {0..ng-1}.
    
    If {ng} is positive and {gv.ne} is nonzero, the pairing is explicit;
    in that case, {gv.ne} must be equal to {ng}, {gv[0]} must be equal
    to {gini}, and {pr[i]} is by definition {gv[i]} for {i} in {0..ng-1}.
*/

/* INTERNAL PROTOTYPES */

bool_t msm_pairing_is_perfect_condensed(msm_pairing_t *pr);
  /* TRUE iff {pr} is a perfect pairing in the condensed format
    (with positive length but empty rung vector). */

void msm_pairing_expand(msm_pairing_t *pr);
  /* If {*pr} is perfect pairing in the condensed format
    (non-empty but with empty rung vector {gv}),
    makes it into an equivalent pairing with explicit
    rungs (non-empty {gv}). */

/* REPRESENTATION-DEPENDENT FUNCTIONS */

int32_t msm_pairing_num_rungs(msm_pairing_t *pr)
  {
    return pr->ng;
  }

bool_t msm_pairing_is_perfect_condensed(msm_pairing_t *pr)
  { return (pr->ng > 0) && (pr->gv.ne == 0); }
  
msm_rung_t msm_pairing_get_rung(msm_pairing_t *pr, int32_t i)
  { /* Get number of rungs {ng}: */
    int32_t ng = pr->ng;
    demand(ng > 0, "empty pairing");
    demand((i >= 0) && (i < ng), "invalid rung index");
    /* Get rung {pr[i]}: */
    msm_rung_t g;
    if (pr->gv.ne == 0)
      { /* Candidate is perfect and condensed, compute rung: */
        g = pr->gini;
        g.c[0] += i; g.c[1] += i;
      }
    else
      { /* Candidate has explicit rungs. */
        assert(pr->gv.ne == ng );
        g = pr->gv.e[i];
      }
    return g;
  }
 
msm_pairing_t *msm_pairing_empty(void)
  { msm_pairing_t *pr = (msm_pairing_t *)notnull(malloc(sizeof(msm_pairing_t)), "no mem");
    (*pr) = (msm_pairing_t)
      { /*ng*/    0,
        /*gini*/  (msm_rung_t){{ 0, 0 }},
        /*gv*/    msm_rung_vec_new(0)
      };
    return pr;
  }
 
msm_pairing_t *msm_pairing_perfect(msm_rung_t gini, int32_t ng)
  { demand(ng > 0, "invalid rung count for perfect pairing");
    /* Compute final indices: */
    msm_pairing_t *pr = (msm_pairing_t *)notnull(malloc(sizeof(msm_pairing_t)), "no mem");
    (*pr) = (msm_pairing_t)
      { /*ng*/    ng,
        /*gini*/  gini,
        /*gv*/    msm_rung_vec_new(0)
      };
    (void)msm_pairing_is_valid(pr, TRUE);
    return pr;
  }

msm_pairing_t *msm_pairing_from_rung_vec(msm_rung_vec_t *gv)
  { int32_t ng = gv->ne;
    demand(ng >= 0, "invalid rung count");
    msm_rung_t gini = (ng > 0 ? gv->e[0] : (msm_rung_t){{ 0, 0 }});
    msm_pairing_t *pr = (msm_pairing_t *)notnull(malloc(sizeof(msm_pairing_t)), "no mem");
    (*pr) = (msm_pairing_t)
      { /*ng*/    ng, 
        /*gini*/  gini,
        /*gv*/    *gv
      };
    (void)msm_pairing_is_valid(pr, TRUE);
    return pr;
  }

void msm_pairing_expand(msm_pairing_t *pr)
  { int32_t ng = pr->ng;
    demand(ng >= 0, "invalid rung count");
    if ((ng > 0) && (pr->gv.ne == 0))
      { /* An implicit perfect pairing: */
        /* Create a new rung vector and fill it with the implicit rungs: */
        msm_rung_vec_t gv = msm_rung_vec_new(ng);
        int32_t i;
        msm_rung_t g = pr->gini;
        for (i = 0; i < ng; i++) { gv.e[i] = g; g.c[0]++; g.c[1]++; }
        pr->gv = gv;
      }
    else
      { /* Already expanded: */
        demand(pr->gv.ne == ng, "inconsistent rung counts");
      }
    /* A bit of paranoia can't hurt: */
    (void)msm_pairing_is_valid(pr, TRUE);
  }

void msm_pairing_free_rungs(msm_pairing_t *pr)
  { if (pr != NULL) 
      { if (pr->gv.e != NULL) 
          { free(pr->gv.e); pr->gv.e = NULL; pr->gv.ne = 0; }
      }
  }

void msm_pairing_free(msm_pairing_t *pr)
  { if (pr != NULL) 
      { if (pr->gv.e != NULL) 
          { free(pr->gv.e); pr->gv.e = NULL; pr->gv.ne = 0; }
        free(pr);
      }
  }
  
void msm_pairing_multi_free(msm_pairing_t *pf[], int32_t maxLevel)  
  { int32_t level;
    for (level = 0; level <= maxLevel; level++)
      { msm_pairing_free(pf[level]); }
  }

/* SEMI-REPRESENTATION-INDEPENDENT FUNCTIONS */
 
msm_pairing_t *msm_pairing_copy(msm_pairing_t *pr)
  { /* Copy the rung vector: */
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    if (msm_pairing_is_perfect_condensed(pr)) 
      { /* Result is perfect condensed too: */
        return msm_pairing_perfect(pr->gini, ng);
      }
    else
      { /* Copy the rungs into a new rung vector: */
        msm_rung_vec_t gv = msm_rung_vec_new(ng);
        int32_t i;
        for (i = 0; i < ng; i++) { gv.e[i] = pr->gv.e[i]; }
        /* Create the pairing: */
        return msm_pairing_from_rung_vec(&gv);
      }
  }
 
msm_pairing_t *msm_pairing_sub_copy(msm_pairing_t *pr, int32_t ini, int32_t fin)
  { demand (ini <= fin, "sub-pairing must be non-empty");
    /* Get the rung count {ng} of {pr}: */
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    
    /* Decide the size {ng_new} of the result's rung vec: */
    int32_t ng_new = fin - ini + 1;
    demand(ini >= 0, "invalid start rung index");
    demand(fin < ng, "invalid end rung index"); 
    
    if (msm_pairing_is_perfect_condensed(pr)) 
      { /* Result is perfect condensed too: */
        msm_rung_t gini_new = msm_pairing_get_rung(pr, ini);
        return msm_pairing_perfect(gini_new, ng_new);
      }
    else
      { /* Copy the relevant rungs into a new rung vector: */
        msm_rung_vec_t gv_new = msm_rung_vec_new(ng_new);
        int32_t i;
        for (i = 0; i < ng_new; i++) { gv_new.e[i] = msm_pairing_get_rung(pr, ini + i); }
        /* Create the pairing: */
        return msm_pairing_from_rung_vec(&(gv_new));
      }
  }

/* REPRESENTATION-INDEPENDENT FUNCTIONS */
 
int32_t msm_pairing_span(msm_pairing_t *pr, int32_t j)
  { int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    if (ng == 0) { return 0; }
    msm_rung_t gini = msm_pairing_get_rung(pr, 0);
    msm_rung_t gfin = msm_pairing_get_rung(pr, ng - 1);
    if (j >= 0)
      { return gfin.c[j] - gini.c[j] + 1; }
    else
      { return 
          (gfin.c[0] - gini.c[0] + 1) +
          (gfin.c[1] - gini.c[1] + 1);
      }
  }
 
msm_pairing_t *msm_pairing_trim_to_span(msm_pairing_t *pr, int32_t j, int32_t span, int32_t i, sign_t dir)
  { demand(span > 0, "cannot trim to zero span");
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng > 0, "cannot trim an empty pairing");
    /* Locate the initial and final rung indices {ini,fin} of the segment: */
    int32_t ini, fin; /* Indices of initial and final rung of segment in {pr}. */
    msm_rung_t gini, gfin; /* Initial and final rung of segment. */
    ini = fin = i;
    gini = gfin = msm_pairing_get_rung(pr, i);
    int32_t segSpan = 1; /* Span of segment between {ini} and {fin} inclusive. */
    bool_t stepped; /* Set to TRUE if the iteration advances. */
    do
      { stepped = FALSE;
        if ((dir <= 0) && (ini > 0))
          { /* Compute span increment {dspan} for decrementing {ini}: */
            msm_rung_t gpre = msm_pairing_get_rung(pr, ini-1);  /* Rung before {ini}. */
            int32_t dspan =  msm_rung_step_span_increment(gpre, gini, j);
            if (segSpan + dspan <= span) 
              { /* Decrement {ini}: */
                ini--; gini = gpre; segSpan += dspan; stepped = TRUE; }
          }
        if ((dir >= 0) && (fin < ng-1))
          { /* Compute span increment {dspan} for decrementing {fin}: */
            msm_rung_t gnex = msm_pairing_get_rung(pr, fin+1); /* Rung after {fin}. */
            int32_t dspan = msm_rung_step_span_increment(gfin, gnex, j);
            if (segSpan + dspan <= span) 
              { /* Increment {fin}: */
                fin++; gfin = gnex; segSpan += dspan; stepped = TRUE;
              }
          }
      }
    while (stepped);
    /* Now extract the segment: */
    return msm_pairing_sub_copy(pr, ini, fin);
  }         
  
int32_t msm_pairing_sub_span(msm_pairing_t *pr, int32_t ini, int32_t fin, int32_t j)
  { if (ini > fin) { return 0; }
    
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng > 0, "pairing is empty");
    demand(ini >= 0, "bad ini");
    demand(fin < ng, "bad fin");
    
    msm_rung_t gini = msm_pairing_get_rung(pr, ini);
    msm_rung_t gfin = msm_pairing_get_rung(pr, fin);
    if (j >= 0)
      { return gfin.c[j] - gini.c[j] + 1; }
    else
      { return 
          (gfin.c[0] - gini.c[0] + 1) +
          (gfin.c[1] - gini.c[1] + 1);
      }
  }

bool_t msm_pairing_equal(msm_pairing_t *pa, msm_pairing_t *pb, bool_t die)
  { int32_t nga = msm_pairing_num_rungs(pa);
    demand(nga >= 0, "invalid {pa} rung count");
    int32_t ngb = msm_pairing_num_rungs(pb);
    demand(ngb >= 0, "invalid {pb} rung count");
    
    if (nga != ngb) 
      { fail_test(die, "rung counts differ"); }
    int32_t t;
    for (t = 0; t < nga; t++)
      { msm_rung_t ga = msm_pairing_get_rung(pa, t);
        msm_rung_t gb = msm_pairing_get_rung(pb, t);
        if (ga.c[0] != gb.c[0]) 
          { fail_test(die, "indices into sequence 0 differ"); }
        if (ga.c[1] != gb.c[1]) 
          { fail_test(die, "indices into sequence 1 differ"); }
      }
    return TRUE;
  }

bool_t msm_pairing_is_increasing(msm_pairing_t *pr, int32_t minIncrEach, int32_t minIncrSum, bool_t atomic, bool_t die)
  { int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    if (ng == 0)
      { /* Empty pairing is increasing and atomic: */
        return TRUE;
      }
    else if (msm_pairing_is_perfect_condensed(pr))
      { /* Condensed perfect pairing is OK if demands are not excessive: */
        return (minIncrEach <= 1) && (minIncrSum <= 2);
      }
    else
      { /* Check the steps in {gv}: */
        msm_rung_t g = msm_pairing_get_rung(pr, 0); /* Previous rung. */
        int32_t k;
        for (k = 1; k < ng ; k++)
          { /* Get next rung: */
            msm_rung_t h = msm_pairing_get_rung(pr, k);
            if (! msm_rung_step_is_increasing(g, h, minIncrEach, minIncrSum, die)) { return FALSE; }
            if (atomic && (! msm_rung_step_is_atomic(g, h, die))) { return FALSE; }
            g = h; 
          }
      }
    return TRUE;
  }
   
bool_t msm_pairing_is_perfect(msm_pairing_t *pr, bool_t die)
  { int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    if (ng == 0)
      { /* Empty pairing, not perfect by definition: */
        return FALSE;
      }
    else if (msm_pairing_is_perfect_condensed(pr))
      { /* Condensed perfect pairing: */
        return TRUE;
      }
    else
      { /* Explicit pairing: must check all steps. */
        msm_rung_t g = pr->gini;
        int32_t k;
        for (k = 1; k < ng ; k++)
          { /* Get next rung: */
            msm_rung_t h = msm_pairing_get_rung(pr, k);
            if (! msm_rung_step_is_perfect(g, h, die)) { return FALSE; }
            g = h; 
          }
        return TRUE;
      }
  }

msm_rung_vec_t msm_pairing_collect_rungs(msm_pairing_t *pr, int32_t n0, int32_t n1)
  { /* Get rung count {ngold} of original pairing: */
    int32_t ng = msm_pairing_num_rungs(pr); /* Number of original rungs. */
    demand(ng >= 0, "invalid rung count");
    /* Reduced rung list is {gv_new[0..ng_new-1]}: */
    msm_rung_vec_t gv_new = msm_rung_vec_new(ng);
    int32_t ng_new = 0; /* Number of distinct reduced rungs. */
    int32_t i;
    for (i = 0; i < ng; i++)
      { /* Get rung number {i} from pairing: */
        msm_rung_t g = msm_pairing_get_rung(pr, i);
        /* Reduce {g} modulo {n0,n1}: */
        if (n0 != 0) { g.c[0] = (int32_t)imod(g.c[0], n0); }
        if (n1 != 0) { g.c[1] = (int32_t)imod(g.c[1], n1); }
        /* Find the proper position {j} of rung {g} in the list {gv_new}: */
        int32_t j = ng_new;
        int32_t cmp = -1;
        while ((j > 0) && ((cmp = msm_rung_lex_compare(gv_new.e[j-1], g)) > 0)) { j--; }
        /* Check wether it is a _new rung: */
        if ((j == ng_new) || (cmp != 0))
          { /* _New rung, insert it in {gv_new}: */
            msm_rung_vec_expand(&gv_new, ng_new);
            if (j < ng_new)
              { /* Open up the slot {gv_new[j]}: */
                int32_t k = ng_new; while (k > j) { gv_new.e[k] = gv_new.e[k-1]; k--; }
              }
            gv_new.e[j] = g;
            ng_new++;
          }
      }
    msm_rung_vec_trim(&gv_new, ng_new);
    return gv_new;
  }

#define msm_rung_huge (msm_rung_t){{ INT32_MAX, INT32_MAX }}
  /* A {msm_rung_t} value that should compare higher than any real rung. */

int32_t msm_pairing_count_common_rungs(msm_pairing_t *pa, msm_pairing_t *pb, int32_t n0, int32_t n1)
  { /* Gather the two rung sets, sorted: */
    msm_rung_vec_t gva = msm_pairing_collect_rungs(pa, n0, n1);
    msm_rung_vec_t gvb = msm_pairing_collect_rungs(pb, n0, n1);
    /* Merge the two lists and count coincidences: */
    int32_t nc = 0;
    int32_t ia = 0, ib = 0;
    while ((ia < gva.ne) || (ib < gvb.ne))
      { msm_rung_t ga = (ia < gva.ne ? gva.e[ia] : msm_rung_huge);
        msm_rung_t gb = (ib < gvb.ne ? gvb.e[ib] : msm_rung_huge);
        int32_t cmp = msm_rung_lex_compare(ga, gb);
        if (cmp < 0)
          { ia++; }
        else if (cmp > 0)
          { ib++; }
        else
          { nc++; ia++; ib++; }
      }
    /* Cleanup: */
    free(gva.e);
    free(gvb.e);
    return nc;
  }

int64_t msm_pairing_perfect_max_length(int32_t n0, int32_t n1)
  { if ((n0 == 0) || (n1 == 0))
      { return 0; }
    else
      { return (n0 < n1 ? n0 : n1); }
  }

void msm_pairing_enum_alignments
  ( int32_t n0,       /* Number of positions on first sequence. */
    int32_t n1,       /* Number of positions on second sequence. */ 
    msm_pairing_perfect_use_proc_t *use
  )
  { fprintf(stderr, "seq 0: %8d positions\n", n0);
    fprintf(stderr, "seq 1: %8d positions\n", n1);
    /* Just in case, check for empty sequences: */
    if ((n0 == 0) || (n1 == 0)) { return; }
    int64_t nrmax = msm_pairing_perfect_max_length(n0, n1); 
    /* Must sweep diagonals of the grid {{0..n0-1}×{0..n1-1}. The diagonals are
      identified by an integer {id} ranging from {idmin} to {idmax}. A
      rung {(i0,i1)} belongs to diagonal number {id = i0 - i1}. */
    int32_t idmin = -(n1 - 1);
    int32_t idmax = +(n0 - 1); 
    fprintf(stderr, "trying %d alignments\n", idmax - idmin + 1);
    int32_t id;
    for(id = idmin; id <= idmax; id++)
      { /* fputc(':', stderr);  */
        /* Computes first index pair {i0,i1} on diagonal: */
        int32_t i0 = (id < 0 ? 0 : +id); 
        int32_t i1 = (id > 0 ? 0 : -id);
        /* Computes number {ng} of index pairs on diagonal: */
        int64_t ng = nrmax;
        if (i0 + ng > n0) { ng = n0 - i0; }
        if (i1 + ng > n1) { ng = n1 - i1; }
        /* Process the alignment: */
        if (ng > 0) { use(i0, i1, ng); }
      }
    fputc('\n', stderr);
  }

int32_t msm_pairing_num_steps(msm_pairing_t *pr)
  {
    return pr->ng-1;
  }

msm_pairing_t *msm_pairing_map_gen(msm_pairing_t *pr, msm_rung_map_proc_t *map)
  { /* Get rung count {ng} of given pairing: */
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    /* Rungs of new pairing will be {gv_new[0..gv_new.ne-1]}: */
    msm_rung_vec_t gv_new = msm_rung_vec_new(ng);
    /* Map all rungs: */
    int32_t i;
    for (i = 0; i < gv_new.ne; i++)
      { gv_new.e[i] = map(msm_pairing_get_rung(pr, i)); }
    /* Make pairing from rung vector: */
    return msm_pairing_from_rung_vec(&gv_new);
  }

msm_pairing_t *msm_pairing_map
  ( msm_pairing_t *pr, 
    msm_seq_desc_t *s0_old,
    msm_seq_desc_t *s1_old,
    msm_seq_desc_t *s0_new,
    msm_seq_desc_t *s1_new
  )
  { 
    auto msm_rung_t map_rung(msm_rung_t g);
    
    msm_rung_t map_rung(msm_rung_t g)
      { double d0 = msm_seq_desc_map_index(g.c[0], s0_old, s0_new);
        double d1 = msm_seq_desc_map_index(g.c[1], s1_old, s1_new);
        int32_t i0 = (int32_t)msm_round(d0);
        int32_t i1 = (int32_t)msm_round(d1);
        return (msm_rung_t){{ i0, i1 }};
      }
      
    return msm_pairing_map_gen(pr, &map_rung);
  }

msm_pairing_t *msm_pairing_interpolate(msm_pairing_t *pr)
  { /* Get number of given rungs {ng}: */
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    /* Get an explicit list of all rungs: */
    msm_rung_vec_t gv = msm_pairing_get_rungs(pr);
    /* Interpolate the rungs: */
    msm_rung_vec_t gv_new = msm_rung_vec_interpolate(&gv);
    /* Reclaim the temporary storage: */
    free(gv.e);
    return msm_pairing_from_rung_vec(&gv_new);
  }
  
msm_pairing_t *msm_pairing_make_increasing(msm_pairing_t *pr, int32_t minIncrEach, int32_t minIncrSum)
  { /* Get number of given rungs {ng}: */
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    /* Get an explicit list of all rungs: */
    msm_rung_vec_t gv = msm_pairing_get_rungs(pr);
    /* Make the rungs increasing: */
    msm_rung_vec_t gv_new = msm_rung_vec_make_increasing(&gv, minIncrEach, minIncrSum);
    /* Reclaim the temporary storage: */
    free(gv.e);
    return msm_pairing_from_rung_vec(&gv_new);
  }

msm_rung_vec_t msm_pairing_get_rungs(msm_pairing_t *pr)
  { /* Get number {ng} of rungs: */
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    /* Allocate the output vector {gv}: */
    msm_rung_vec_t gv = msm_rung_vec_new(ng);
    /* Copy all rungs to {gv}: */
    int32_t i;
    for (i = 0; i < ng; i++) { gv.e[i] = msm_pairing_get_rung(pr, i); }
    return gv;
  }

msm_pairing_t *msm_pairing_throw
  ( int32_t n0,
    int32_t n1, 
    int32_t len,
    bool_t atomic,
    bool_t diagonal,
    double skipProb
  )
  { demand((n0 > 0) && (n1 > 0), "cannot pair up an empty sequence");
    demand(len > 0, "cannot pair up an empty sequence");

    msm_pairing_t *pr = (msm_pairing_t *)notnull(malloc(sizeof(msm_pairing_t)), "no mem");
    /* Shall we generate a perfect pairing? */
    bool_t perfect = ((atomic && (skipProb == 0.0)) || (drandom() < 0.1));
    if (perfect)
      { /* Generate an implicit perfect pairing. */
        /* Compute the max number of rungs {ngmax}: */
        int64_t ngmax = msm_pairing_perfect_max_length(n0, n1);
        
        int32_t ng = (int32_t)imin(ngmax, len);
        /* Compute the maximum initial rung {gmax}: */
        msm_rung_t gmin = (msm_rung_t){{ 0, 0 }};
        msm_rung_t gmax = (msm_rung_t){{ n0 - ng, n1 - ng }};
        /* Choose the initial rung {gini}: */
        msm_rung_t gini = msm_rung_throw(gmin, gmax, diagonal);
        /* Generate pairing: */
        pr = msm_pairing_perfect(gini, ng);
      }
    else
      { /* Generate an explicit pairing. */
        /* Choose the initial rung {gini} and final rung {gfin}. */
        msm_rung_t gini; /* Initial rung. */
        msm_rung_t gfin; /* Final rung. */
        /* Choose a segment of sequence 0: */
        int32_t ini0, fin0;
        msm_seq_desc_throw_segment(n0, len, &ini0, &fin0);
        /* Coose a segment of sequence 1: */
        int32_t ini1, fin1;
        msm_seq_desc_throw_segment(n1, len, &ini1, &fin1);
        /* The pairing should span both segments: */
        gini = (msm_rung_t){{ ini0, ini1 }};
        gfin = (msm_rung_t){{ fin0, fin1 }};
        
        /* Now generate a random path from {gini} to {gfin}: */
        msm_rung_vec_t gv = msm_rung_vec_throw(gini, gfin, atomic, skipProb);
       /* Wrap it up as a pairing: */
        pr = msm_pairing_from_rung_vec(&gv);
      }
    return pr;
  }

msm_pairing_t *msm_pairing_sub_throw(msm_pairing_t *pr, int32_t len)
  {
    int32_t ng = msm_pairing_num_rungs(pr);
    int32_t span = 2*len; /* Required total span. */
    
    /* Find the index {kMax} in {0..ng-1} of the last rung that can be initial: */
    int32_t kMax;
    {int32_t k = 0;
    while ((k < ng) && (msm_pairing_sub_span(pr, k, ng-1, -1) >= span)) { k++; }
    kMax = k - 1;}
    
    /* Pick the index {ini} in {0..kMax} of the central rung of the pairing: */
    int32_t ini = int32_abrandom(0, kMax);
    
    /* Find the index {fin} of the last rung of the pairing: */
    int32_t fin;
    { int32_t k = ini + 1;
      while (msm_pairing_sub_span(pr, ini, k, -1) < span) { k++; }
      fin = k - 1;
    }
    
    /* Extract the rungs between {ini} and {fin}, inclusive: */
    msm_rung_vec_t gv = msm_rung_vec_new(fin - ini + 1);
    int32_t k;
    for (k = ini; k <= fin; k++)
      { gv.e[k - ini] = msm_pairing_get_rung(pr, k); }
    return msm_pairing_from_rung_vec(&gv);
  }

void msm_pairing_write_full(FILE *wr, msm_pairing_t *pr)
  {
    int32_t ng = msm_pairing_num_rungs(pr);
    int32_t k;
    fprintf(wr, "%d rungs\n", ng);
    for (k = 0; k < ng; k++)
      { msm_rung_t g = msm_pairing_get_rung(pr, k);
        fprintf(wr, "  %08d ( %08d , %08d )\n", k, g.c[0], g.c[1]);
      }
  }

void msm_pairing_write(FILE *wr, msm_pairing_t *pr, int32_t ixSize, bool_t writePairing)
  {
    /* Write number of pairs {ng}: */
    int32_t ng = msm_pairing_num_rungs(pr);
    demand(ng >= 0, "invalid rung count");
    fprintf(wr, "%*d", ixSize, ng);
    
    /* Write initial rung: */
    msm_rung_t gini = (ng == 0 ? (msm_rung_t){{ 0, 0 }} : msm_pairing_get_rung(pr, 0));
    fprintf(wr, " %*d %*d", ixSize, gini.c[0], ixSize, gini.c[1]);

    if(writePairing){
      fprintf(wr, " ");
      /* Print rung vector: */
      if ((ng == 0) || msm_pairing_is_perfect(pr, FALSE))
        { /* Perfect pairing (condensed or not) or empty pairing; omit the rung vector: */ 
          fprintf(wr, "*");
        }
      else
        { /* Imperfect pairing;  write the perfect sub-pairings explicitly: */ 
          /* Pretend the step before the first rung was perfect so that the first rung is '|': */ 
          msm_rung_t g = (msm_rung_t){{ gini.c[0] - 1, gini.c[1] - 1 }}; /* Previous rung. */
          msm_rung_t h = gini; /* Next rung to analyze, or {msm_rung_none}. */
          int32_t i = 0; 
          while (i < ng)
            { /* Get the gaps {d0,d1} from the previous rung: */
              int32_t d0 = h.c[0] - g.c[0];
              int32_t d1 = h.c[1] - g.c[1];
              /* Count the rungs {rep} until the next non-perfect step of end of pairing: */
              int32_t rep = 0;
              do
                { g = h;
                  i++; rep++;
                  if (i >= ng) { h = msm_rung_none; break; }
                  h = msm_pairing_get_rung(pr, i);
                }
              while ((h.c[0] - g.c[0] == 1) && (h.c[1] - g.c[1] == 1));
              /* Write the rung set: */
              msm_rung_step_write(wr, d0, d1, rep);
            }
        }
    }
    fflush(wr);
  }

msm_pairing_t *msm_pairing_read(FILE *rd)
  {
    bool_t debug = FALSE;
    
    /* Number of rungs: */
    int32_t ng = fget_int32(rd);
    demand(ng >= 0, "bad number of pairs");
    
    /* Initial rung: */
    msm_rung_t gini;
    gini.c[0] = fget_int32(rd); demand(gini.c[0] >= 0, "bad init index for sequence 0"); 
    gini.c[1] = fget_int32(rd); demand(gini.c[1] >= 0, "bad init index for sequence 1");
    if (ng == 0)
      { demand((gini.c[0] == 0) && (gini.c[1] == 0), "empty pairing with nonzero initial rung"); }
    
    /* Read first char of pairing: */
    fget_skip_spaces(rd); 
    int32_t ch = fget_char(rd);
    
    if (ch == '*')
      { if (ng == 0)
          { /* Empty pairing: */
            return msm_pairing_empty();
          }
        else
          { /* Non-empty perfect pairing: */
            return msm_pairing_perfect(gini, ng);
          }
      }
    else if (ch == '|')
      { /* Generic pairing. */
        demand(ng > 0, "empty pairing with explicit rungs");
        /* Allocate the rung vector: */
        msm_rung_vec_t gv = msm_rung_vec_new(ng);
        /* Fake previous rung: */
        msm_rung_t g = (msm_rung_t){{ gini.c[0] - 1, gini.c[1] - 1 }};
        /* Read the perfect sub-pairings: */
        int32_t kg = 0; /* Next free position in {gv}. */
        while (TRUE)
          { demand (ch != EOF, "unexpected EOF");
            ungetc(ch, rd);
            if 
              ( (ch != '\'') && (ch != '.') && (ch != ':') && (ch != '(') && 
                (ch != '|') && (ch != '\\') && (ch != '/') 
              ) { break; }
            /* Read another perfect sub-pairing: */
            int32_t d0, d1, rep;
            msm_rung_step_read(rd, &d0, &d1, &rep);
            if (debug) { fprintf(stderr, "    read: d0 = %8d  d1 = %8d  rep = %8d", d0, d1, rep); }
            assert((d0 >= 0) && (d1 >= 0));
            assert(rep > 0);
            if (kg == 0) { assert((d0 == 1) && (d1 == 1)); }
            /* Append it: */
            int32_t r;
            for (r = 0; r < rep; r++) 
              { demand(kg < ng, "too many pairs in pairing desc");
                g.c[0] += d0; d0 = 1;
                g.c[1] += d1; d1 = 1;
                gv.e[kg] = g; 
                kg++;
              }
            if (debug) { fprintf(stderr, "  kg = %d\n", kg); }
            /* Skip to next perfect sub-pairing: */
            fget_skip_spaces(rd); 
            ch = fgetc(rd);
          }
        /* Now {rd} is positioned just before the next non-pairing char: */
        assert(kg <= ng);
        demand(kg == ng, "too few pairs in pairing desc");
        return msm_pairing_from_rung_vec(&gv);
      }
    else
      { demand(FALSE, "bad or missing initial pairing character");
        return NULL; /* Compiler pacifier. */
      }
  }

bool_t msm_pairing_is_valid (msm_pairing_t *pr, bool_t die)
  { int32_t ng = pr->ng;
    if (ng < 0) { fail_test(die, "invalid rung count"); }
    /* Validity of {gv}: */
    if (pr->gv.ne == 0) 
      { /* Empty pairing or condensed perfect pairing, {ng} may be zero or more: */
        if (pr->gv.e != NULL) { fail_test(die, "non-null vector with zero length"); }
      }
    else
      { /* Explicit (non-condensed) pairing: */
        if (pr->gv.ne != ng) 
          { fail_test(die, "inconsistent rung vector length"); }
        if (pr->gv.e == NULL) 
          { fail_test(die, "null address of non-empty rung vector"); }
        /* First rung must be consistent: */
        if (pr->gini.c[0] != pr->gv.e[0].c[0]) 
          { fail_test(die, "first rung is inconsistent on side 0"); }
        if (pr->gini.c[1] != pr->gv.e[0].c[1]) 
          { fail_test(die, "first rung is inconsistent on side 1"); }
        /* Steps must be non-decreasing: */
        msm_rung_t g = pr->gini;
        int32_t i;
        for (i = 1; i < ng; i++) 
          { msm_rung_t h = pr->gv.e[i];
            if (! msm_rung_step_is_increasing(g, h, 0, 0, die)) { return FALSE; }
            g = h;
          }
      }
    return TRUE;
  }
 
