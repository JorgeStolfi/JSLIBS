/* See msm_pairing.h */
/* Last edited on 2008-01-12 08:26:48 by stolfi */

#define msm_pairing_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <assert.h>
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <sign.h>
#include <fget.h>
#include <nget.h>
#include <jsrandom.h>
#include <jsmath.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>
#include <msm_rung.h>
#include <msm_pairing.h>

/* INTERNAL PAIRING REPRESENTATION */

struct msm_pairing_t 
  { int ng;              /* Number of pairs in pairing, or rung index period. */
    int circp;           /* TRUE iff pairing is circular. */
    msm_rung_t gini;      /* Rung number 0. */
    msm_rung_vec_t gv;    /* Rungs of pairing, or empty if perfect. */
  };
  /* The rung count {ng} is always positive.
    
    If {gv.nel} is 0, the pairing is implicit and perfect -- namely,
    {p[i] = gini + (i,i)} for all {i} in {0..ng-1}. If {gv.nel != 0},
    the pairing is explicit; in that case, {gv.nel} must be equal to
    {ng}, {gv[0]} must be equal to {gini}, and {p[i]} is by
    definition {gv[i]} for {i} in {0..ng-1}.

    If {circp} is false, the pairing is open and contains the rungs
    {p[0..ng-1]} only. If {circp} is true, the pairing is circular; in
    that case, the fundamental cycle is {p[0..ng-1]}, and, for all
    {i}, rung {p[i+ng]} is {p[i] + D}, where {D = p[ng] - p[0]}.
    
    If {circp} is true and the pairing is not implicit, then {gv.nel}
    must be {ng+1}, and {p[ng]} is by definition {gv.el[ng]}. */

/* INTERNAL PROTOTYPES */

bool_t msm_pairing_is_perfect_condensed(msm_pairing_t *p);
  /* TRUE iff {p} is a perfect pairing in the condensed format
    (with empty rung vector). */

void msm_pairing_expand(msm_pairing_t *p);
  /* If {*p} is an implicit perfect pairing (with empty rung vector {gv}),
    makes it into an equivalent generic pairing (with non-empty {gv}). */

/* REPRESENTATION-DEPENDENT FUNCTIONS */

vec_typeimpl(msm_rung_vec_t,msm_rung_vec,msm_rung_t);

bool_t msm_pairing_is_perfect_condensed(msm_pairing_t *p)
  { return p->gv.nel == 0; }

int msm_pairing_num_rungs(msm_pairing_t *p)
  { return p->ng; }
  
bool_t msm_pairing_is_circular(msm_pairing_t *p)
  { return p->circp; }

msm_rung_t msm_pairing_period(msm_pairing_t *p)
  { int ng = p->ng;
    if (p->circp)
      { if (p->gv.nel == 0)
          { return (msm_rung_t){{ ng, ng }}; }
        else
          { assert(p->gv.nel == ng+1);
            msm_rung_t grep = p->gv.el[ng];
            return (msm_rung_t){{ grep.c[0] - p->gini.c[0],  grep.c[1] - p->gini.c[1] }};
          }
      }
    else
      { return (msm_rung_t){{ 0, 0 }}; }
  }
  
msm_rung_t msm_pairing_get_rung(msm_pairing_t *p, int i)
  { /* Get number of rungs {ng}: */
    int ng = p->ng;
    demand(ng > 0, "invalid rung count");
    /* Get circularity flag {circp}: */
    bool_t circp = msm_pairing_is_circular(p);
    /* Reduce the index {i} to a rung index {k} in the fundamental period: */ 
    int k = (circp ? imod(i, ng) : i); 
    demand((k >= 0) && (k < ng), "invalid rung index");
    /* Get rung {p[k]}: */
    msm_rung_t g;
    if (p->gv.nel == 0)
      { /* Candidate is perfect and condensed, compute rung: */
        g = p->gini;
        g.c[0] += k; g.c[1] += k;
      }
    else
      { /* Candidate has explicit rungs. */
        assert(p->gv.nel == ng + (circp ? 1 : 0));
        g = p->gv.el[k];
      }
    /* If circular, add period vectors as needed: */
    if (circp)
      { int q = (i - k)/ng; /* Number of periods between {k} and {i}. */
        msm_rung_t per = msm_pairing_period(p);
        g.c[0] += q*per.c[0]; g.c[1] += q*per.c[1]; 
      }
    return g;
  }

bool_t msm_pairing_is_valid (msm_pairing_t *p, bool_t die)
  { int ng = p->ng;
    if (ng <= 0) { fail_test(die, "invalid rung count"); }
    /* Validity of {gv}: */
    if (p->gv.nel < 0) { fail_test(die, "invalid rung vector length"); }
    if (p->gv.nel == 0) 
      { /* Implicit perfect pairing: */
        if (p->gv.el != NULL) { fail_test(die, "non-null vector with zero length"); }
      }
    else
      { /* Explicit pairing: */
        if (p->gv.nel != ng + (p->circp ? 1 : 0)) 
          { fail_test(die, "inconsistent rung vector length"); }
        if (p->gv.el == NULL) 
          { fail_test(die, "null address of non-empty rung vector"); }
        /* First rung must be consistent: */
        if (p->gini.c[0] != p->gv.el[0].c[0]) 
          { fail_test(die, "first rung is inconsistent in X"); }
        if (p->gini.c[1] != p->gv.el[0].c[1]) 
          { fail_test(die, "first rung is inconsistent in X"); }
      }
    return TRUE;
  }
  
msm_pairing_t *msm_pairing_perfect(msm_rung_t *gini, int ng, bool_t circp)
  { /* A pairing must have at least one rung: */
    demand(ng > 0, "invalid rung count for perfect pairing");
    /* Compute final indices: */
    msm_pairing_t *p = (msm_pairing_t *)notnull(malloc(sizeof(msm_pairing_t)), "no mem");
    (*p) = (msm_pairing_t)
      { /*ng*/    ng,
        /*circp*/ circp,
        /*gini*/  *gini,
        /*gv*/    msm_rung_vec_new(0)
      };
    (void)msm_pairing_is_valid(p, TRUE);
    return p;
  }

msm_pairing_t *msm_pairing_from_rung_vec(msm_rung_vec_t *gv, bool_t circp)
  { int ng = gv->nel - (circp ? 1 : 0);
    demand(ng > 0, "invalid rung count");
    msm_pairing_t *p = (msm_pairing_t *)notnull(malloc(sizeof(msm_pairing_t)), "no mem");
    (*p) = (msm_pairing_t)
      { /*ng*/    ng, 
        /*circp*/ circp,
        /*gini*/  gv->el[0],
        /*gv*/    *gv
      };
    (void)msm_pairing_is_valid(p, TRUE);
    return p;
  }

void msm_pairing_expand(msm_pairing_t *p)
  { int ng = p->ng;
    if (p->gv.nel == 0)
      { /* Must be an implicit perfect pairing: */
        assert(ng > 1);
        /* Create a new rung vector and fill it with a perfect pairing: */
        p->gv = msm_rung_vec_new(ng);
        int i;
        msm_rung_t g = p->gini;
        for (i = 0; i < ng; i++)
          { p->gv.el[i] = g; g.c[0]++; g.c[1]++; }
      }
    /* A bit of paranoia can't hurt: */
    (void)msm_pairing_is_valid(p, TRUE);
  }

void msm_pairing_free(msm_pairing_t *p)
  { if (p != NULL) 
      { if (p->gv.el != NULL) 
          { free(p->gv.el); p->gv.el = NULL; p->gv.nel = 0; }
        free(p);
      }
  }

/* REPRESENTATION-INDEPENDENT FUNCTIONS */
  
bool_t msm_pairing_equal(msm_pairing_t *pa, msm_pairing_t *pb, bool_t die)
  { int nga = msm_pairing_num_rungs(pa);
    int ngb = msm_pairing_num_rungs(pb);
    if (nga != ngb) 
      { fail_test(die, "rung counts differ"); }
    bool_t circpa = msm_pairing_is_circular(pa);
    bool_t circpb = msm_pairing_is_circular(pb);
    if (circpa != circpb) 
      { fail_test(die, "circularities differ"); }
    int t;
    for (t = 0; t < nga; t++)
      { msm_rung_t ga = msm_pairing_get_rung(pa, t);
        msm_rung_t gb = msm_pairing_get_rung(pb, t);
        if (ga.c[0] != gb.c[0]) 
          { fail_test(die, "X coordinates of rungs differ"); }
        if (ga.c[1] != gb.c[1]) 
          { fail_test(die, "Y coordinates of rungs differ"); }
      }
    return TRUE;
  }

bool_t msm_pairing_is_increasing(msm_pairing_t *p, bool_t die)
  { int ng = msm_pairing_num_rungs(p);
    demand(ng > 0, "invalid rung count");
    if (msm_pairing_is_perfect_condensed(p))
      { /* Condensed perfect pairing: */
        return TRUE;
      }
    else
      { /* Check the fundamental steps in {gv}: */
        bool_t circp = msm_pairing_is_circular(p);
        msm_rung_t h = msm_pairing_get_rung(p, 0);
        int k;
        for (k = 1; k < ng + (circp ? 1 : 0); k++)
          { msm_rung_t g = h; 
            /* Get next rung: */
            h = msm_pairing_get_rung(p, k);
            if (! msm_rung_step_is_increasing(&g, &h, die)) { return FALSE; }
          }
      }
    return TRUE;
  }

bool_t msm_pairing_is_atomic_increasing(msm_pairing_t *p, bool_t die)
  { int ng = msm_pairing_num_rungs(p);
    demand(ng > 0, "invalid rung count");
    if (msm_pairing_is_perfect_condensed(p))
      { /* Condensed perfect pairing: */
        return TRUE;
      }
    else
      { /* Check steps in {gv}: */
        bool_t circp = msm_pairing_is_circular(p);
        msm_rung_t h = msm_pairing_get_rung(p, 0);
        int k;
        for (k = 1; k < ng + (circp ? 1 : 0); k++)
          { msm_rung_t g = h; 
            /* Get next rung (or first one if circular pairing): */
            h = msm_pairing_get_rung(p, k);
            if (! msm_rung_step_is_increasing(&g, &h, die)) { return FALSE; }
            if (! msm_rung_step_is_atomic(&g, &h, die)) { return FALSE; }
          }
      }
    return TRUE;
  }
   
bool_t msm_pairing_is_perfect(msm_pairing_t *p, bool_t die)
  { int ng = msm_pairing_num_rungs(p);
    demand(ng > 0, "invalid rung count");
    if (msm_pairing_is_perfect_condensed(p))
      { /* Condensed perfect pairing: */
        return TRUE;
      }
    else
      { /* Explicit pairing: must check all steps. */
        bool_t circp = msm_pairing_is_circular(p);
        msm_rung_t g = p->gini;
        int k;
        for (k = 1; k < ng + (circp ? 1 : 0); k++)
          { /* Get next rung: */
            msm_rung_t h = msm_pairing_get_rung(p, k);
            if (! msm_rung_step_is_perfect(&g, &h, die)) { return FALSE; }
            g = h; 
          }
        return TRUE;
      }
  }

msm_rung_vec_t msm_pairing_collect_rungs(msm_pairing_t *p, int nx, int ny)
  { /* The coordinate periods of a circular pairing must be multiples of {nx,ny}: */
    if (msm_pairing_is_circular(p))
      { msm_rung_t per = msm_pairing_period(p);
        if (nx != 0) { demand(per.c[0] % nx == 0, "X period inconsistent"); }
        if (ny != 0) { demand(per.c[1] % ny == 0, "Y period inconsistent"); }
      }
    /* Get rung count {ngold} of original pairing: */
    int ngold = msm_pairing_num_rungs(p); /* Number of original rungs. */
    /* Reduced rung list is {rgvnew[0..ngnew-1]}: */
    msm_rung_vec_t rgvnew = msm_rung_vec_new(ngold);
    int ngnew = 0; /* Number of distinct reduced rungs. */
    int i;
    for (i = 0; i < ngold; i++)
      { /* Get rung number {i} from pairing: */
        msm_rung_t g = msm_pairing_get_rung(p, i);
        /* Reduce {g} modulo {nx,ny}: */
        if (nx != 0) { g.c[0] = imod(g.c[0], nx); }
        if (ny != 0) { g.c[1] = imod(g.c[1], ny); }
        /* Find the proper position {j} of rung {g} in the list {rgvnew}: */
        int j = ngnew;
        int cmp;
        while ((j > 0) && ((cmp = msm_rung_lex_compare(&(rgvnew.el[j-1]), &g)) > 0)) { j--; }
        /* Check wether it is a new rung: */
        if ((j == ngnew) || (cmp != 0))
          { /* New rung, insert it in {rgvnew}: */
            msm_rung_vec_expand(&rgvnew, ngnew);
            if (j < ngnew)
              { /* Open up the slot {rgvnew[j]}: */
                int k = ngnew; while (k > j) { rgvnew.el[k] = rgvnew.el[k-1]; k--; }
              }
            rgvnew.el[j] = g;
          }
      }
    return rgvnew;
  }

#define msm_rung_huge (msm_rung_t){{ INT_MAX, INT_MAX }}
  /* A {msm_rung_t} value that should compare higher than any real rung. */

int msm_pairing_count_common_rungs(msm_pairing_t *pa, msm_pairing_t *pb, int nx, int ny)
  { /* Gather the two rung sets, sorted: */
    msm_rung_vec_t gva = msm_pairing_collect_rungs(pa, nx, ny);
    msm_rung_vec_t gvb = msm_pairing_collect_rungs(pb, nx, ny);
    /* Merge the two lists and count coincidences: */
    int nc = 0;
    int ia = 0, ib = 0;
    while ((ia < gva.nel) || (ib < gvb.nel))
      { msm_rung_t ga = (ia < gva.nel ? gva.el[ia] : msm_rung_huge);
        msm_rung_t gb = (ib < gvb.nel ? gvb.el[ib] : msm_rung_huge);
        int cmp = msm_rung_lex_compare(&ga, &gb);
        if (cmp < 0)
          { ia++; }
        else if (cmp > 0)
          { ib++; }
        else
          { nc++; ia++; ib++; }
      }
    /* Cleanup: */
    free(gva.el);
    free(gvb.el);
    return nc;
  }

int64_t msm_pairing_perfect_max_length(int nx, bool_t circx, int ny, bool_t circy)
  { if ((nx == 0) || (ny == 0))
      { return 0; }
    else if (circx && circy)
      { return lcm(nx, ny); }
    else if (circx)
      { return ny; }
    else if (circy)
      { return nx; }
    else
      { return (nx < ny ? nx : ny); }
  }

void msm_pairing_enum_alignments
  ( int nx,       /* Number of bases on first sequence. */
    bool_t circx, /* TRUE if first sequence is circular. */
    int ny,       /* Number of bases on second sequence. */ 
    bool_t circy, /* TRUE if second sequence is circular. */ 
    msm_pairing_perfect_use_proc_t *use
  )
  { fprintf(stderr, "x seq: %8d bases circp = %d\n", nx, circx);
    fprintf(stderr, "y seq: %8d bases circp = %d\n", ny, circy);
    bool_t circp = (circx & circy); /* TRUE iff pairings are circular. */
    /* Just in case, check for empty sequences: */
    if ((nx == 0) || (ny == 0)) { return; }
    int idmin, idmax;
    int64_t nrmax = msm_pairing_perfect_max_length(nx, circx, ny, circy); 
    int id;
    /* Must sweep diagonals of the {xp × yp} grid. The diagonals are
      identified by an integer {id} ranging from {idmin} to {idmax}.
      A pair of samples with indices {ix} and {iy} belongs to diagonal
      number {id = ix - iy}. Note that if either sequence is circular
      the diagonals will wrap around. */
    if (circx && circy)
      { idmin = 0; idmax = gcd(nx, ny) - 1; }
    else if (circx)
      { idmin = 0; idmax = +(nx - 1); }
    else if (circy)
      { idmin = -(ny - 1); idmax = 0; }
    else
      { idmin = -(ny - 1); idmax = +(nx + 1); }
    fprintf(stderr, "trying %d alignments\n", idmax - idmin + 1);
    for(id = idmin; id <= idmax; id++)
      { /* fputc(':', stderr);  */
        /* Computes range of pair index {jd} along diagonal: */
        int ix, iy;
        int64_t ng;
        if (circx && circy)
          { ix = +id; iy = 0; ng = nrmax; }
        else if (circx)
          { ix = +id; iy = 0; ng = nrmax; }
        else if (circy)
          { ix = 0; iy = -id; ng = nrmax; }
        else
          { ix = (id < 0 ? 0 : +id); 
            iy = (id > 0 ? 0 : -id);
            ng = nrmax;
            if (ix + ng > nx) { ng = nx - ix; }
            if (iy + ng > ny) { ng = ny - iy; }
          }    
        /* Process the alignment: */
        if (ng > 0) { use(ix, iy, ng, circp); }
      }
    fputc('\n', stderr);
  }

msm_pairing_t *msm_pairing_map_coordinates(msm_pairing_t *p, msm_rung_map_proc_t *map)
  { /* Get rung count {ng} and circularity {circp} of given pairing: */
    int ng = msm_pairing_num_rungs(p);
    bool_t circp = msm_pairing_is_circular(p);
    /* Rungs of new pairing will be {gvnew[0..gvnew.nel-1]}: */
    msm_rung_vec_t gvnew = msm_rung_vec_new(ng + (circp ? 1 : 0));
    /* Map all rungs (including the loopback, if circular): */
    int i;
    for (i = 0; i < gvnew.nel; i++)
      { gvnew.el[i] = map(msm_pairing_get_rung(p, i)); }
    /* Make pairing from rung vector: */
    return msm_pairing_from_rung_vec(&gvnew, circp);
  }

msm_pairing_t *msm_pairing_interpolate(msm_pairing_t *p)
  { /* Get number of given rungs {ng} and circularity: */
    int ng = msm_pairing_num_rungs(p);
    assert(ng > 0);
    bool_t circp = msm_pairing_is_circular(p);
    /* Get an explicit list of all rungs (including the loopback rung if circular): */
    msm_rung_vec_t gvold = msm_pairing_get_rungs(p);
    /* Interpolate the rungs: */
    msm_rung_vec_t gvnew = msm_rung_vec_interpolate(&gvold);
    /* Reclaim the temporary storage: */
    free(gvold.el);
    return msm_pairing_from_rung_vec(&gvnew, circp);
  }
  
msm_rung_vec_t msm_pairing_get_rungs(msm_pairing_t *p)
  { /* Get number {ng} of fundamental rungs: */
    int ng = msm_pairing_num_rungs(p);
    assert(ng > 0);
    bool_t circp = msm_pairing_is_circular(p);
    /* Get number {ngout} of output rungs and allocate the output vector {gvout}: */
    int ngout = (circp ? ng+1 : ng);
    msm_rung_vec_t gv = msm_rung_vec_new(ngout);
    /* Copy all rungs to {gvout}: */
    int i;
    for (i = 0; i < ngout; i++) 
      { gv.el[i] = msm_pairing_get_rung(p, i);
        /* Paranoid check that step {g->h} is strictly increasing: */
        if (i > 0) { (void)msm_rung_step_is_increasing(&(gv.el[i-1]), &(gv.el[i]), TRUE); }
      }
    return gv;
  }

msm_pairing_t *msm_pairing_throw
  ( int nx,
    bool_t circx, 
    int ny, 
    bool_t circy, 
    int len,
    bool_t circp, 
    bool_t atomic,
    bool_t diagonal,
    double skipProb
  )
  { demand((nx > 0) && (ny > 0), "cannot pair up an empty sequence");
    if (circp) { demand(circx && circy, "a circular pairing requires circular seqs"); }
    msm_pairing_t *p = (msm_pairing_t *)notnull(malloc(sizeof(msm_pairing_t)), "no mem");
    /* Shall we generate a perfect pairing? */
    bool_t perfect = ((atomic && (skipProb == 0.0)) || (drandom() < 0.1));
    if (perfect)
      { /* Generate an implicit perfect pairing. */
        /* Compute the max number of rungs {ngmax}: */
        int ngmax = msm_pairing_perfect_max_length(nx, circx, ny, circy);
        /* Choose the number of rungs {ng}.
          A circular perfect pairing must have precisely {ngmax} rungs.
          An open perfect pairing may have less than that. */
        int ng = (circp ? ngmax : len);
        /* Compute the maximum initial rung {gmax}: */
        msm_rung_t gmin = (msm_rung_t){{ 0, 0 }};
        msm_rung_t gmax = (msm_rung_t){{ nx - (circx ? 1 : ng), ny - (circy ? 1 : ng) }};
        /* Choose the initial rung {gini}: */
        msm_rung_t gini = msm_rung_throw(&gmin, &gmax, diagonal);
        /* Generate pairing: */
        p = msm_pairing_perfect(&gini, ng, circp);
      }
    else
      { /* Generate an explicit pairing. */
        /* Choose the initial rung {gini} and final (or loopback) rung {gfin}. */
        msm_rung_t gini; /* Initial rung. */
        msm_rung_t gfin; /* Final rung (loopback if {circp}). */
        if (circp)
          { msm_rung_t gmin = (msm_rung_t){{ 0, 0 }};
            msm_rung_t gmax = (msm_rung_t){{ nx-1, ny-1 }};
            gini = msm_rung_throw(&gmin, &gmax, diagonal);
            /* The coord periods must be multiples of the seq lengths.
              In most applications it makes no sense to have pairings 
              that make funny numbers of turns around each seq.
              So we generate a 1-increasing pairing that does 
              {ny} turns around the X-sequence and {nx} turns around
              the Y-sequence:
            */
            int DXY = lcm(nx, ny);
            gfin = (msm_rung_t){{ gini.c[0] + DXY, gini.c[1] + DXY }};
          }
        else
          { /* The pairing is open. */
            /* Coose a segment of the X sequence (possibly with wrap-around): */
            int xini, xfin;
            msm_seq_desc_throw_segment(nx, circx, len, &xini, &xfin);
            /* Coose a segment of the Y sequence (possibly with wrap-around): */
            int yini, yfin;
            msm_seq_desc_throw_segment(ny, circy, len, &yini, &yfin);
            /* The pairing should span both segments: */
            gini = (msm_rung_t){{ xini, yini }};
            gfin = (msm_rung_t){{ xfin, yfin }};
          }
        /* Now generate a random path from {gini} to {gfin}: */
        msm_rung_vec_t gv = msm_rung_vec_throw(&gini, &gfin, atomic, skipProb);
       /* Wrap it up as a pairing: */
        p = msm_pairing_from_rung_vec(&gv, circp);
      }
    return p;
  }

void msm_pairing_write(FILE *wr, msm_pairing_t *p, int ixsize)
  {
    /* Write number of pairs {ng}: */
    int ng = msm_pairing_num_rungs(p);
    demand(ng > 0, "invalid {ng}");
    fprintf(wr, "%*d", ixsize, ng);
    
    /* Write circularity flag: */
    bool_t circp = msm_pairing_is_circular(p);
    fprintf(wr, " %c", (circp ? 'T' : 'F'));

    /* Write initial rung: */
    msm_rung_t gini = msm_pairing_get_rung(p, 0);
    fprintf(wr, " %*d %*d", ixsize, gini.c[0], ixsize, gini.c[1]);

    /* Print rung vector: */
    if (msm_pairing_is_perfect(p, FALSE))
      { /* Perfect pairing (condensed or not); omit the rung vector: */ 
        fprintf(wr, " *");
      }
    else
      { /* Imperfect pairing;  write the rungs explicitly: */ 
        /* Print code of first pair (always "|"). */
        fprintf(wr, " |");
        /* Print codes of other pairs: */
        msm_rung_t g = gini; /* Last printed rung. */
        int i; 
        for (i = 1; i < ng + (circp ? 1 : 0); i++)
          { msm_rung_t h = msm_pairing_get_rung(p, i);
            msm_rung_step_write(wr, &g, &h);
            g = h;
          }
      }
    fflush(wr);
  }

msm_pairing_t *msm_pairing_read(FILE *rd)
  {
    /* Number of rungs: */
    int ng = fget_int(rd);
    demand(ng > 0, "bad number of pairs");
    
    /* Circularity flag: */
    bool_t circp = fget_bool(rd);
    
    /* Initial rung: */
    msm_rung_t gini;
    gini.c[0] = fget_int(rd); demand(gini.c[0] >= 0, "bad init index for {x} sequence"); 
    gini.c[1] = fget_int(rd); demand(gini.c[1] >= 0, "bad init index for {y} sequence"); 
    
    /* Read first char of pairing: */
    fget_skip_spaces(rd); 
    int ch = fget_char(rd);
    
    if (ch == '*')
      { /* Perfect matching: */
        return msm_pairing_perfect(&gini, ng, circp);
      }
    else if (ch == '|')
      { /* Generic pairing. */
        /* Allocate the rung vector: */
        int gvsize = ng + (circp ? 1 : 0);
        msm_rung_vec_t gv = msm_rung_vec_new(gvsize);
        /* Save initial pair: */
        msm_rung_t g = gini;
        gv.el[0] = g;
        /* Read and save remaining pairs: */
        int i; 
        for(i = 1; i < gv.nel; i++)
          { msm_rung_t h = msm_rung_step_read(rd, g); 
            gv.el[i] = h;
            g = h;
          }
        return msm_pairing_from_rung_vec(&gv, circp);
      }
    else
      { demand(FALSE, "bad or missing initial pairing character");
        return NULL; /* Compiler pacifier. */
      }
  }

msm_pairing_t *msm_pairing_map_to_finer
  ( msm_pairing_t *p, 
    int nxold, 
    int nyold, 
    int nxnew, 
    int nynew, 
    int nwtb
  )
  { 
    auto msm_rung_t map_rung(msm_rung_t g);
    
    msm_rung_t map_rung(msm_rung_t g)
      { int ix = msm_seq_desc_map_index_to_finer(g.c[0], nxold, nxnew, nwtb);
        int iy = msm_seq_desc_map_index_to_finer(g.c[1], nyold, nynew, nwtb);
        return (msm_rung_t){{ix,iy}};
      }
      
    return msm_pairing_map_coordinates(p, &map_rung);
  }

msm_pairing_t *msm_pairing_map_to_coarser
  ( msm_pairing_t *p, 
    int nxold, 
    int nyold, 
    int nxnew, 
    int nynew, 
    int nwtb
  )
  { 
    auto msm_rung_t map_rung(msm_rung_t g);
    
    msm_rung_t map_rung(msm_rung_t g)
      { int ix = msm_seq_desc_map_index_to_coarser(g.c[0], nxold, nxnew, nwtb);
        int iy = msm_seq_desc_map_index_to_coarser(g.c[1], nyold, nynew, nwtb);
        return (msm_rung_t){{ix,iy}};
      }
      
    return msm_pairing_map_coordinates(p, &map_rung);
  }

