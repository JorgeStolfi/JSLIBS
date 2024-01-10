/* See msm_cand_vec.h */
/* Last edited on 2011-12-22 01:16:21 by stolfi */

#define msm_cand_vec_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <jswsize.h>

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>
#include <msm_cand_refine.h>
#include <msm_dyn.h>

#define msm_cand_vec_checking_level 1
  /* Define this as positive to get paranoid checking. */

#define CHKLEV(n) (msm_cand_vec_checking_level >= (n))
  /* Use {if(CHKLEV(n)){...}} for checks of paranoia level {n}. */

#define msm_cand_vec_debug_level 2
  /* Set this positive to get some debugging printouts. */

#define DEBLEV(n) (msm_cand_vec_debug_level >= (n))
  /* Use {if(DEBLEV(n)){...}} for debug printouts of level {n}. */

vec_typeimpl(msm_cand_vec_t,msm_cand_vec,msm_cand_t);

typedef double msm_subpairing_score_proc_t(int64_t r, int64_t s, bool_t circsub); 
  /* A client procedure that computes a score for a subpairing {p[r..s]}
    of some perfect pairing {p}.  If {circsub} the sub-pairing is circular
    and equal to the whole of {p},  which should be circular too. */

typedef void msm_subpairing_use_proc_t
  ( int64_t r,
    int64_t s,
    bool_t circsub,
    double score
  ); 
  /* A client procedure that uses a subpairing {p[r..s]} of some
    perfect pairing {p}, whose score is {score}. If {circsub} the
    subpairing is circular and equal to the whole of {p}, which should
    be circular too. */

void msm_cand_scan_alignment
  ( int64_t ng,         /* Number of pairs in alignment. */
    bool_t circp,        /* TRUE if alignment is circular. */ 
    msm_subpairing_score_proc_t *score,
    int64_t minRungs,   /* Min number of rungs in a valid subpairing. */
    int64_t maxRungs,   /* Max number of rungs in a valid subpairing. */
    msm_subpairing_use_proc_t *use  /* What to do with candidates. */
  );
  /* Scans the perfect pairing {p} of two sequences {x,y} that has
    {ng} rungs, and calls {use} on all candidates found in it. If
    {circp} the pairing {p} is assumed to be circular.

    A valid sub-pairing is a set of consecutive rungs from {p} that
    has at least {minRungs} and at most {maxRungs} rungs. A candidate
    is a valid sub-pairing whose score is positive and better than
    that of any of proper valid sub-pairing or super-pairing.
    
    The score of a non-circular subpairing {p[r..s]} is assumed to be
    given by {score(r,s,FALSE)}. If {p} is circular, the score of {p}
    should be given by {score(0,ng-1,TRUE)}. */

msm_cand_vec_t msm_cand_vec_get_best_perfect
  ( msm_seq_desc_t *xp,
    msm_seq_desc_t *yp,
    msm_rung_step_score_proc_t *step_score,
    int64_t minRungs, 
    double minScore,
    int maxCands,
    double frac
  )
  {
    int nx = xp->npos; bool_t circx = xp->circ;
    int ny = yp->npos; bool_t circy = yp->circ;

    int64_t maxRungs = (nx < ny ? nx : ny);  
      /* Max number of rungs in a candidate, including the circular one. */
    
    /* Working vectors for {msm_cand_scan_alignment}: */
    double_vec_t rsc = double_vec_new(0);
    double_vec_t psc = double_vec_new(0);

    msm_cand_vec_t cdv = msm_cand_vec_new(maxCands);
    int ncd = 0;
    /* The current set of candidates is kept in the array
      {cdv[0..ncd-1]}. During the scan, the array is sorted
      by *decreasing* score (i.e. the worst candidate is 
      {cdv[ncd-1]}. */

    auto void scan_alignment(int ix, int iy, int64_t ng, bool_t circp);
      /* Scans the perfect pairing {p} that starts with rung {(ix,iy)}
        and has {ng} rungs, and adds all candidates found in it to the
        heap. If {circp} the pairing is assumed to be circular.
        
        A valid sub-pairing is a set of consecutive rungs from {p}
        that has at least {minRungs} and at most {min(nx,ny)} rungs. A
        candidate is a valid sub-pairing whose score is positive and
        better than that of any of proper valid sub-pairing or
        super-pairing. */

    /* Scan all maximal-length perfect pairings, collect best cands: */
    msm_pairing_enum_alignments(nx, circx, ny, circy, &scan_alignment); 
    
    /* Release working storage: */
    double_vec_trim(&rsc, 0);
    double_vec_trim(&psc, 0);

    /* Trim the candidate vector to the number of cands actually found: */
    msm_cand_vec_trim(&cdv, ncd);
    
    return cdv; 

    void scan_alignment(int ix, int iy, int64_t ng, bool_t circp)
      {
        demand(ng > 0, "bad {ng}");

        auto double subpairing_score(int64_t r, int64_t s, bool_t circsub);
          /* Score of sub-pairing {p[r..s]}, namely
            {S = SUM{step_score(p[j-1],p[j]) : j \in r+1..s}}.
            If {circsub} adds also {step_score(p[s],p[r])},
            else adds also none-to-first-rung and last-rung-to-none. 
            Obtained from precomputed tables.
          */

        auto double pref_score(int64_t k);
          /* Score of non-circular sub-pairing {p[0..k]}, from tables. */

        auto double join_score(int64_t k);
          /* Score of rung {p[k]}, from tables. */

        auto void initialize_score_tables(void);  
          /* Initializes the tables used by {score}. */
          
        auto void check_all_scores(void); 
          /* Checks {subpairing_score(r,s)} against defintion for all valid {(r,s)}. */
         
        auto void check_score(int64_t r, int64_t s);  
          /* Checks {subpairing_score(r,s)} against defintion. */

        auto void store_candidate(int64_t r, int64_t s, bool_t circsub, double score);
          /* Stores into the sorted list {cdv} the pairing
            {(ix+k,iy+k)} for {k=r,..s}, where both indices are taken
            modulo the respective sequence lengths. If {circsub},
            assumes that the pairing is circular. Assumes that its
            score is {score}. The procedure is a no-op if the pairing
            is not among the best {maxCands} entries found so far, or
            is already in the list. */
        
        if (DEBLEV(1))
          { fprintf(stderr, ("\nalignment (%d,%d:%" int64_d_fmt ")\n"), ix, iy, ng); }
            
        /* Allocate score tables: */
        initialize_score_tables();

        /* Paranoia checks for the {subpairing_score} function: */
        if (CHKLEV(2)) { check_all_scores(); }

        msm_cand_scan_alignment
          ( ng, circp, 
            &subpairing_score,
            minRungs, 
            maxRungs, 
            &store_candidate
          );
          
        auto void update_score_tables(int64_t kcur);
          /* Update the table of precomputed scores so that it is valid for {k} 
            up to {kcur}.  May discard some older entries. */

        int64_t kmin, kmax, ksize;
        /* The tables {rsc,psc} are used to store precomputed rung and
          total scores for this alignment. The score of rung {p[k]} is
          stored in {rsc[imod(k,ksize)]}. The score of the (non-circular)
          sub-pairing {p[0..k]} is stored in {psc[imod(k,ksize)]}. These data
          are valid only for {k} in {kmin..kmax}. */

        void initialize_score_tables(void)
          { ksize = 2*maxRungs - minRungs + 1;
            double_vec_expand(&rsc, ksize);
            double_vec_expand(&psc, ksize);
            msm_rung_t g = (msm_rung_t){{ix, iy}};
            rsc.e[0] = step_score(xp, yp, msm_rung_none, g);
            psc.e[0] = rsc.e[0];
            kmin = 0; kmax = 0; 
          }

        void update_score_tables(int64_t kcur)
          { 
            demand(kcur >= kmin, "table is too short");
            if (kcur <= kmax) { return; }
            msm_rung_t g = (msm_rung_t){{ix + kmax, iy + kmax}};
            do
              { int64_t k = kmax + 1;
                msm_rung_t h = (msm_rung_t){{ix + k, iy + k}};
                double ssck = (step_score == NULL ? 0 : step_score(xp, yp, g, h));
                int kk = imod(k,ksize);
                psc.e[kk] = psc.e[imod(kmax,ksize)] + ssck;
                g = h;
                kmax++;
                if ((kmax - kmin) > ksize) { kmin++; }
              }
            while (kmax < kcur);
          }

        double pref_score(int64_t k)
          { demand(k >= 0, "invalid prefix: k < 0");
            demand(circp || (k < ng), "invalid prefix: out of range");
            /* Table consistency: */
            update_score_tables(k);
            affirm((kmin <= k) && (k <= kmax), "out of table range");
            return psc.e[imod(k,ksize)]; 
          }

        double join_score(int64_t k)
          { demand(k >= 0, "invalid prefix: k < 0");
            demand(circp || (k < ng), "invalid prefix: out of range");
            /* Table consistency: */
            update_score_tables(k);
            affirm((kmin <= k) && (k <= kmax), "out of table range");
            return rsc.e[imod(k,ksize)]; 
          }

        double subpairing_score(int64_t a, int64_t b, bool_t circsub)
          { /* Argument consistency: */
            /* fprintf(stderr, ("  sc a = %6" int64_d_fmt "  b = %6" int64_d_fmt " %c"), a, b, (circsub ? 'C' : 'N')); */
            demand(a <= b, "invalid sub-pairing: a > b");
            demand(a >= 0, "invalid sub-pairing: a < 0");
            demand(circp || (b < ng), "invalid sub-pairing: out of range");
            demand(b - a + 1 <= maxRungs, "invalid sub-pairing: too long");
            double s;
            if (circsub)
              { assert((b - a) == ng - 1);
                s = pref_score(b + 1) - pref_score(a);
              }
            else
              { s = pref_score(b) - pref_score(a) + join_score(a); }
            /* fprintf(stderr, " = %24.16e\n", s); */
            return s;
          }
          
        void check_score(int64_t r, int64_t s)
          { double eps = 1.0e-10; /* Absolute error allowed. */
            double tol = 1.0e-6;  /* Relative error allowed. */
            msm_rung_t gini = (msm_rung_t){{ix+r, iy+r}};
            msm_pairing_t *pf = msm_pairing_perfect(gini, s-r+1, FALSE);
            double s_slow = msm_cand_compute_score(xp, yp, pf, step_score);
            double s_fast = subpairing_score(r,s, FALSE);
            double s_size = fabs(s_slow) + fabs(s_fast) + eps;
            if (fabs(s_slow - s_fast) > eps + tol*s_size)
              { fprintf(stderr, "bad {score}\n");
                fprintf(stderr, "  fast = %24.26e\n", s_fast);
                fprintf(stderr, "  slow = %24.26e\n", s_slow);
                assert(FALSE);
              }
            msm_pairing_free(pf);
          }

        void check_all_scores(void)
          { fprintf(stderr, "checking the {subpairing_score} function...\n");
            int64_t rt, st;
            for (rt = 0; rt < ng; rt++)
              { 
                int64_t stmax = rt + maxRungs - 1;
                if ((! circp) && (stmax >= ng)) { stmax = ng - 1; }
                for (st = rt + minRungs - 1; st <= stmax; st++)
                  {  check_score(rt, st); }
                fprintf(stderr, "#"); 
              }
            fprintf(stderr, "\n"); 
          }

        void store_candidate(int64_t r, int64_t s, bool_t circsub, double score)
          { /* Reduce the initial indices modulo the respective periods: */
            msm_rung_t gini = (msm_rung_t){{ imod(ix + r, nx), imod(iy + r, ny) }};
            int64_t ns = s - r + 1;
            if (DEBLEV(2))
              { fprintf
                  ( stderr, 
                    ("\n  perfect candidate starting at (%d,%d) steps = %" int64_d_fmt " score = %24.16e\n"), 
                    gini.c[0], gini.c[1], ns, score
                  ); 
              }
            /* Worth adding: */
            if (score < minScore) { return; }
            /* Compute the coordinate period vector: */
            if (circsub) { demand(ns == lcm(nx,ny), "invalid perfect circular pairing"); }
            msm_pairing_t *pr = msm_pairing_perfect(gini, ns, circsub);
            msm_cand_t cd = msm_cand_from_pairing(xp, yp, pr, score);
            msm_cand_vec_insert(&cdv, &ncd, 1, maxCands, &cd, frac);
          }
      }
  }

void msm_cand_scan_alignment
  ( int64_t ng,         /* Number of pairs in alignment. */
    bool_t circp,       /* TRUE if alignment is circular. */ 
    msm_subpairing_score_proc_t *score,
    int64_t minRungs,   /* Min number of rungs in a valid subpairing. */
    int64_t maxRungs,   /* Max number of rungs in a valid subpairing. */
    msm_subpairing_use_proc_t *use  /* What to do with candidates. */
  )
  {
    /* Trivial cases: */
    if ((ng < minRungs) || (maxRungs < minRungs)) { return; }
    
    auto void initialize_mxsma(int64_t ri, int64_t si);
      /* Initialize the {mxsma} table for the current subpairing {[ri,si]}: */
    
    auto void initialize_mxbig(int64_t ri, int64_t si);
      /* Initialize the {mxbig} table for the current subpairing {[ri,si]}: */

    auto double *mxsma(int64_t r, int64_t s);
      /* The maximum score among all valid sub-pairings {[a,b]} that 
        are contained in or equal to {[r,s]}. */
      
    auto double *mxbig(int64_t r, int64_t s);
      /* The maximum score among all valid  subpairings {[a,b]} 
        that properly contain {[r,s]}. */
    
    auto void update_mxsma(int64_t rn, int64_t sn);
      /* Recompute {mxsma(t,sn)} from {mxsma(t,sn-1)} 
        for {t} in {rn .. sn-minRungs+1}. */

    auto void update_mxbig(int64_t rn, int64_t sn);
      /* Recompute {mxbig(rn,t)} from {mxbig(rn-1,t)} 
        for {t} in {sn .. rn + maxRungs - 1}. */

    auto void check_mxsma(int64_t rn, int64_t sn);
      /* Checks {mxsma(t,sn)} for {t} in {rn .. sn-minRungs+1}. */

    auto void check_mxbig(int64_t rn, int64_t sn);
      /* Checks {mxbig(rn,t)} for {t} in {sn .. rn + maxRungs - 1}. */

    int N = maxRungs - minRungs + 1;
    double tsma[N]; /* The value of {mxsma(r,s)} is stored in {tsma[imod(r,N)]} */
    double tbig[N]; /* The value of {mxbig(r,s)} is stored in {tbig[imod(s,N)]} */

    /* If the pairing {p} is circular, we must consider it too as a
      potential candidate. In that case, note that {p} properly
      contains any other candidate; so, if its score is better than
      that of any open sub-pairing, {p} is the only candidate. In all
      other cases {p} is not a candidate. */
      
    bool_t circ_valid = circp && (ng <= maxRungs); 
      /* TRUE if the full circular pairing is valid. */
    
    double circ_score = (circ_valid ? score(0, ng, TRUE) : -INF);
      /* Score of full circular pairing, if it is valid, else {-oo}. */
    
    bool_t circ_best = circ_valid;
      /* TRUE while the full curcular alignment is the best one. */

    /* Find the initial pre-candidate {[r,s]}.
      
      Pre-candidate condition: A valid sub-pairing {[r,s]} is a
      pre-candidate if {s} is the smallest index such that {[r,s]} is
      valid and {mxsma(r,s)> mxbig(r,s)}. A pre-candidate is a
      candidate if {f(r,s)==mxsma(r,s)}.  */
    
    int64_t rc, sc; /* Current pre-candidate is {[rc,sc]}. */
    
    /* Choose {rc} of initial pre-candidate: */
    if (! circp)
      { /* Must start at {rc = 0}. Note that {mxbig} will work fine. */
        rc = 0;
      }
    else
      { /* Start at a large enough {sc} so that {mxbig(rc,sc)} 
          can be computed without going into negative indices. */
        rc = maxRungs - minRungs;
      }
      
    /* Intialize {sc} at minimum value for this {rc}: */
    sc = rc + minRungs - 1;
      
    /* Initialize the {mxsma,mxbig} tables for this subpairing: */
    initialize_mxsma(rc,sc);
    initialize_mxbig(rc,sc);
    
    /* Scan pre-candidates and output candidates: */
    int64_t ncol;
    for (ncol = 0; ncol < ng - 1; ncol++)
      { 
        /* fprintf(stderr, ("\nlooking for candidates with r = %" int64_d_fmt "\n"), rc); */
        /* fprintf(stderr, ("  starting with [%" int64_d_fmt ",%" int64_d_fmt "]\n"), rc, sc); */
        /* Get max {sc} for this {rc}: */
        int64_t smax = rc + maxRungs - 1;
        if ((! circp) && (smax >= ng)) { smax = ng - 1; }
        assert(sc <= smax);
        /* Increment {sc} trying to satisfy the pre-candidate condition: */
        double vsma = (*(mxsma(rc,sc)));
        double vbig = (*(mxbig(rc,sc)));
        /* fprintf(stderr, "  mxsma = %24.16e  mxbig = %24.16e\n", vsma, vbig); */
        while ((sc < smax) && (vsma <= vbig)) 
          { sc++;
            /* fprintf(stderr, ("  incrementing sc to [%" int64_d_fmt ",%" int64_d_fmt "]\n"), rc, sc); */
            update_mxsma(rc,sc);
            vsma = (*(mxsma(rc,sc)));
            vbig = (*(mxbig(rc,sc)));
            /* fprintf(stderr, ("  mxsma = %24.16e  mxbig = %24.16e\n"), vsma, vbig); */
          }
        /* If {rc,sc} is a candidate, output it: */
        if (vsma > vbig)
          { double rssc = score(rc,sc,FALSE);
            if ((rssc > circ_score) && (rssc == vsma)) 
              { use(rc, sc, FALSE, rssc);
                circ_best = FALSE; /* The full circular alignment is not the best one. */
              }
          }
        /* Advance {rc}, incrementing {sc} if necessary to keep {rc,sc} valid: */
        int64_t ns = sc - rc + 1;
        if (ns <= minRungs)
          { assert(ns == minRungs);
            rc++; 
            sc++; 
            if ((! circp) && (sc >= ng)) { /* No more cands. */  break; }
            /* fprintf(stderr, ("  incrementing rc and sc to [%" int64_d_fmt ",%" int64_d_fmt "]\n"), rc, sc); */
            update_mxsma(rc,sc);
            update_mxbig(rc,sc);
          }
        else
          { rc++;
            /* fprintf(stderr, ("  incrementing rc to [%" int64_d_fmt ",%" int64_d_fmt "]\n"), rc, sc); */
            update_mxbig(rc,sc);
          }
      }
    
    /* If full circular subpairing is still the best, add it: */
    if (circ_best) { use(0, ng-1, TRUE, circ_score); }

    void initialize_mxsma(int64_t ri, int64_t si)
      { int64_t rmax = si - minRungs + 1;
        assert(rmax >= ri);
        /* Compute {mxsma(t,si)} for all {t} ii {ri..rmax}: */
        /* fprintf(stderr, ("  initializing mxsma(%" int64_d_fmt "..%" int64_d_fmt ",%" int64_d_fmt ")\n"), ri, rmax, si); */
        int64_t smin = ri + minRungs - 1;
        int64_t s1;
        for (s1 = smin; s1 <= si; s1++) { update_mxsma(ri, s1); }
      }
          
    void initialize_mxbig(int64_t ri, int64_t si)
      {
        int64_t smax = ri + maxRungs - 1;
        if ((! circp) && (smax >= ng)) { smax = ng; }
        assert(smax >= si);
        /* Compute {mxbig(ri,t)} for {t} ii {si..smax}: */
        /* fprintf(stderr, ("  initializing mxbig(%" int64_d_fmt ",%" int64_d_fmt "..%" int64_d_fmt ")\n"), ri, si, smax); */
        int rmin = si - maxRungs + 1;
        if ((! circp) && (rmin < 0)) { rmin = 0; }
        assert(rmin >= 0); /* If {circ} caller must use {ri} large enough to avoid negs. */
        int64_t r1;
        for (r1 = rmin; r1 <= ri; r1++) { update_mxbig(r1, si); }
      }

    double *mxsma(int64_t r, int64_t s)
      { assert(r >= 0); return &(tsma[r % N]); }
      
    double *mxbig(int64_t r, int64_t s)
      { assert(s >= 0); return &(tbig[s % N]); }
      
    void update_mxsma(int64_t rn, int64_t sn)
      { assert(rn >= 0);
        assert(sn >= 0);
        double scmax = -INF;
        int64_t rmax = sn - minRungs + 1; 
        assert(rmax >= rn);
        int64_t rt;
        for (rt = rmax; rt >= rn; rt--)
          { 
            double s;
            s = score(rt, sn, FALSE);
            if (s > scmax) { scmax = s; }
            if (rt < rmax)
              { s = (*(mxsma(rt,sn-1)));
                if (s > scmax) { scmax = s; }
              }
            /* fprintf(stderr, ("  updating mxsma(%" int64_d_fmt ",%" int64_d_fmt ") = %24.16e\n"), rt, sn, scmax); */
            (*(mxsma(rt,sn))) = scmax;
          }
        if (CHKLEV(1)) { check_mxsma(rn, sn); }
      }
      
    void check_mxsma(int64_t rn, int64_t sn)
      { assert(rn >= 0);
        assert(sn >= 0);
        double scmax = -INF;
        int64_t rmax = sn - minRungs + 1; 
        assert(rmax >= rn);
        /* fprintf(stderr, ("  checking mxsma(%" int64_d_fmt "..%" int64_d_fmt ",%" int64_d_fmt ")\n"), rn, rmax, sn); */
        int64_t rt;
        for (rt = rmax; rt >= rn; rt--)
          { int64_t smin = rt + minRungs - 1;
            assert(smin >= 0);
            assert(smin <= sn);
            int64_t st;
            for (st = smin; st <= sn; st++)
              { double s = score(rt, st, FALSE);
                if (s > scmax) { scmax = s; }
              }
            assert((*(mxsma(rt,sn))) == scmax);
          }
      }
      
    void update_mxbig(int64_t rn, int64_t sn)
      { assert(rn >= 0);
        assert(sn >= 0);
        double scmax = -INF;
        int64_t smax = rn + maxRungs - 1; 
        if ((! circp) && (smax >= ng)) { smax = ng - 1; }
        assert(smax >= sn);
        int64_t st;
        for (st = smax; st >= sn; st--)
          { 
            int64_t rmin = st - maxRungs + 1;
            if ((! circp) && (rmin < 0)) { rmin = 0; }
            assert(rmin >= 0); /* If {p} is circular, {rn} must be large enough to avoid negs. */
            double s;
            if (rn > rmin)
              { s = score(rn-1, st, FALSE);
                if (s > scmax) { scmax = s; }
                s = (*(mxbig(rn-1,st)));
                if (s > scmax) { scmax = s; }
              }
            /* fprintf(stderr, ("  updating mxbig(%" int64_d_fmt ",%" int64_d_fmt ") = %24.16e\n"), rn, st, scmax); */
            (*(mxbig(rn,st))) = scmax;
            s = score(rn, st, FALSE);
            if (s > scmax) { scmax = s; }
          }
        if (CHKLEV(1)) { check_mxbig(rn, sn); }
      }

    void check_mxbig(int64_t rn, int64_t sn)
      { assert(rn >= 0);
        assert(sn >= 0);
        double scmax = -INF;
        int64_t smax = rn + maxRungs - 1; 
        if ((! circp) && (smax >= ng)) { smax = ng - 1; }
        assert(smax >= sn);
        /* fprintf(stderr, ("  checking mxbig(%" int64_d_fmt ",%" int64_d_fmt "..%" int64_d_fmt ")\n"), rn, sn, smax); */
        int64_t st;
        for (st = smax; st >= sn; st--)
          { int64_t rmin = st - maxRungs + 1;
            if ((! circp) && (rmin < 0)) { rmin = 0; }
            assert(rmin >= 0); /* If {p} is circular, {rn} must be large enough to avoid negs. */
            int64_t rt;
            double s;
            for (rt = rmin; rt < rn; rt++)
              { s = score(rt, st, FALSE);
                if (s > scmax) { scmax = s; }
              }
            assert((*(mxbig(rn,st))) == scmax);  
            s = score(rn, st, FALSE);
            if (s > scmax) { scmax = s; }
          }
      }
  }

msm_cand_vec_t msm_cand_vec_throw
  ( int ncd,
    msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp,
    int minlen,
    int maxlen,
    double circProb,
    double atomProb,
    double diagProb,
    double skipProb
  )
  { /* Allocate the candidate vector {cdv}: */
    msm_cand_vec_t cdv = msm_cand_vec_new(ncd);
    int i; 
    for (i = 0; i < ncd; i++) 
      { /* Choose between circular and open pairing: */
        bool_t circP = xp->circ & yp->circ & (drandom() < circProb);
        /* Choose between 1-atomic or not: */
        bool_t atomP = (drandom() < atomProb);
        /* Choose between near-diagonal and random: */
        bool_t diagP = (drandom() < diagProb);
        /* Choose length of candidate: */
        int len = abrandom(minlen, maxlen);
        /* Generate random candidate: */
        cdv.e[i] = msm_cand_throw(xp, yp, len, circP, atomP, diagP, 0.10);
      }
    return cdv;
  }

void msm_cand_vec_free(msm_cand_vec_t *cdv)
  {
    int i;
    for (i = 0; i < cdv->ne; i++) { msm_cand_free(&(cdv->e[i])); }
    free(cdv->e);
  }

#define msm_cand_vec_type_name "msm_cand_vec"
#define msm_cand_vec_format "2008-01-11"

void msm_cand_vec_write(FILE *wr, msm_cand_vec_t *cdv)
  { int ncd = cdv->ne;
    /* Find maximum abs of {id,strlen(name),index}: */
    int idmax = 1; 
    int nlenmax = 1;
    int ixmax = 1;
    int i;
    for (i = 0; i < ncd; i++)
      { msm_cand_t *cd = &(cdv->e[i]);
        int j;
        for (j = 0; j < 2; j++)
          { msm_seq_desc_t *sj = &(cd->seq[j]);
            int idj = (sj->id >= 0 ? sj->id : 10*(-sj->id)); 
            if (idj > idmax) { idmax = idj; }
            int lenj = (sj->name == NULL ? 0 : strlen(sj->name));
            if (lenj > nlenmax) { nlenmax = lenj; }
            int ixj = (sj->npos <= 0 ? 0 : sj->npos - 1);
            if (ixj > ixmax) { ixmax = ixj; }
          }
      }
    /* Compute the max number of digits needed for {ix,name,index}: */
    int idSize = digits(idmax);
    int nameSize = nlenmax;
    int ixSize = digits(ixmax);
    filefmt_write_header(wr, msm_cand_vec_type_name, msm_cand_vec_format);
    fprintf(wr, "ncands = %d\n", ncd);
    for (i = 0; i < ncd; i++)
      { msm_cand_write(wr, NULL, &(cdv->e[i]), idSize, nameSize, ixSize, NULL);
        fputc('\n', wr);
      }
    filefmt_write_footer(wr, msm_cand_vec_type_name);
    fflush(wr);
  }

void msm_cand_vec_write_named(msm_cand_vec_t *cdv, char *name, char *tag)
  { FILE *wr = msm_open_write(name, tag, ".cdv", TRUE);
    msm_cand_vec_write(wr, cdv);
    fclose(wr);
  }
  
msm_cand_vec_t msm_cand_vec_read(FILE *rd)
  { /* Check and skip the file header: */
    filefmt_read_header(rd, msm_cand_vec_type_name, msm_cand_vec_format);
    /* Skip comment lines, if any: */
    (void)filefmt_read_comment(rd, '|');
    /* Read the header fields: */
    bool_t ncd = nget_int(rd, "ncands"); fget_eol(rd);
    /* Read the candidates: */
    msm_cand_vec_t cdv = msm_cand_vec_new(ncd);
    int k; 
    for (k = 0; k < ncd; k++)
      { cdv.e[k] = msm_cand_read(rd, NULL, NULL); 
        fget_eol(rd);
      }
    /* Check and skip file footer: */
    filefmt_read_footer(rd, msm_cand_vec_type_name);
    return cdv;
  }
 
msm_cand_vec_t msm_cand_vec_read_named(char *name, char *tag)
  { FILE *rd = msm_open_read(name, tag, ".cdv", TRUE);
    msm_cand_vec_t cdv = msm_cand_vec_read(rd);
    fclose(rd);
    return cdv;
  }

msm_cand_vec_t msm_cand_vec_map_to_finer
  ( msm_cand_vec_t *cdv,
    int nxnew,
    int nynew,
    bool_t circx,
    bool_t circy,
    int nwtb
  )
  { fprintf(stderr, "  mapping %d candidates ...\n", cdv->ne);
    msm_cand_vec_t cdvnew = msm_cand_vec_new(cdv->ne);
    int ic;
    for (ic = 0; ic < cdv->ne; ic++)
      { cdvnew.e[ic] = msm_cand_map_to_finer
          ( &(cdv->e[ic]), nxnew, nynew, circx, circy, nwtb );
      }
    fprintf(stderr, "  mapped %d candidates.\n", cdvnew.ne);
    return cdvnew;
  }
 
msm_cand_vec_t msm_cand_vec_interpolate(msm_cand_vec_t *cdv)
  { msm_cand_vec_t cdvnew = msm_cand_vec_new(cdv->ne);
    int ic;
    for (ic = 0; ic < cdv->ne; ic++)
      { cdvnew.e[ic] = msm_cand_interpolate(&(cdv->e[ic])); }
    return cdvnew;
  }

msm_cand_vec_t msm_cand_vec_refine
  ( msm_cand_vec_t *cdvold,
    msm_seq_desc_t *ap, 
    msm_seq_desc_t *bp,
    int delta,
    int kappa, 
    bool_t expand,
    bool_t shrink,
    int maxUnp, 
    msm_rung_step_score_proc_t *step_score,
    msm_dyn_tableau_t *tb, 
    int minCover,
    int maxCands,
    double frac
  )
  { msm_cand_vec_t cdvnew = msm_cand_vec_new(cdvold->ne);
    int ic;
    int ncnew = 0; /* Number of candidates that were kept. */
    for (ic = 0; ic < cdvold->ne; ic++)
      { msm_cand_t *cdoi = &(cdvold->e[ic]);
        demand(msm_seq_desc_equal(&(cdoi->seq[0]), ap, FALSE), "wrong seq 0");
        demand(msm_seq_desc_equal(&(cdoi->seq[1]), bp, FALSE), "wrong seq 1");
        msm_cand_t cdni = msm_cand_refine
          ( cdoi, delta, kappa, expand, shrink, maxUnp, step_score, tb );
        msm_cand_vec_insert(&cdvnew, &ncnew, minCover, maxCands, &cdni, frac);
      }
    msm_cand_vec_trim(&cdvnew, ncnew);
    return cdvnew;
  }

void msm_cand_vec_insert
  ( msm_cand_vec_t *cdv,
    int *ncandP,
    int minCover,
    int maxCands,
    msm_cand_t *cd,
    double frac
  )
  {
    auto int bubble_up(int j);
      /* Assumes that the list is OK except that {cdv[j]} may be worse
        than {cdv[j-1]}. Swaps the two. */
    
    auto void delete(int j);
      /* Deletes the entry {cdv[j]}, displacing all subsequent entries. 
        Decrements {ncand}. */
    
    /* Check minimum coverage requirement: */
    int j;
    for (j = 0; j < 2; j++)
      { if (msm_pairing_span(cd->pr, j) < minCover) 
        { /* Candidate does not cover enough of sequence {j}: */
          msm_cand_debug("N", cd);
          return;
        }
      }
    
    int k; /* Final position of inserted candidate, or {(*ncandP)} if not inserted. */

    /* Find the smallest {k} such that {cdv[k]} is definitely worse than {cd}: */
    k = (*ncandP);
    while ((k > 0) && (cdv->e[k-1].score < cd->score)) { k--; }

    if (frac >= 0.0)
      { /* Discard {cd} if there is any better or equal cand that contains or is similar to {cd}: */
        for (j = 0; j < k; j++)
          { int ngnew = msm_pairing_fund_rung_count(cd->pr); /* Rungs in {cd}. */
            int ngboth = msm_cand_count_common_rungs(&(cdv->e[j]), cd); /* Shared rungs. */
            if (ngboth >= frac*ngnew) 
              { /* Discard {cd}: */
                msm_cand_debug("C", cd);
                return;
              }
          }

        /* Remove any candidates that are worse than {cd} but contained in or similar to {cd}: */
        for (j = (*ncandP)-1; j >= k; j--)
          { int ngold = msm_pairing_fund_rung_count(cdv->e[j].pr); /* Rungs in {cdv->e[j]}. */
            int ngboth = msm_cand_count_common_rungs(&(cdv->e[j]), cd); /* Shared rungs. */
            if (ngboth >= frac*ngold) 
              { /* Discard candidate {cvd.e[j]}: */
                msm_cand_debug("S", &(cdv->e[j]));
                delete(j);
              }
          }
      }

    if ((*ncandP) >= maxCands)
      { k = (*ncandP)-1;
        if (cd->score < cdv->e[k].score)
          { /* Discard {cd}: */
            msm_cand_debug("B", cd);
            return;
          }
        else
          { /* Delete worst candidate: */
            msm_cand_debug("-", &(cdv->e[k]));
            delete(k);
          }
      }

    /* Add {cd} to {cdv[0..ncand-1]}: */
    assert((*ncandP) < maxCands);
    k = (*ncandP);
    msm_cand_vec_expand(cdv, k);
    cdv->e[k] = (*cd);
    (*ncandP)++;
    k = bubble_up(k);
    assert((k >= 0) && (k < (*ncandP)));
    msm_cand_debug("+", &(cdv->e[k]));

    int bubble_up(int j)
      { while (j > 0)
          { /* Get parent: */
            int k = j-1; 
            /* If {cdv[k]} is not worse than {cdv[j]}, stop: */
            if (cdv->e[k].score >= cdv->e[j].score) { break; }
            /* Swap parent with {cdv[j]}: */
            msm_cand_t pp = cdv->e[k]; cdv->e[k] = cdv->e[j]; cdv->e[j] = pp;
            j = k;
          }
        return j;
      }

    void delete(int j)
      { int i;
        for (i = j+1; i < (*ncandP); i++) { cdv->e[i-1] = cdv->e[i]; }
        (*ncandP)--;
      }
  }
