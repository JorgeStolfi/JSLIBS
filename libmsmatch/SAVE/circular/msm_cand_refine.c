/* See msm_cand_refine.h */
/* Last edited on 2008-04-19 16:22:22 by stolfi */

#define msm_refine_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <assert.h>
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <filefmt.h>
#include <vec.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>
#include <msm_rung.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_cand_refine.h>

/* INTERNAL PROTOTYPES */

void msm_cand_refine_fill_tableau
  ( msm_cand_t *cdold, 
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp,            /* Maximum unpaired datums between rungs. */
    msm_rung_step_score_proc_t *step_score, 
    msm_rung_t *goptp,     /* OUT: last rung of optimum pairing. */
    double *voptp,         /* OUT: score of optimum pairing. */
    msm_dyn_tableau_t *tb  /* WORK: dynamic programming tableau. */
  );
  /* Runs the incremental dynamic programming algorithm to compute an
    optimum pairing of two sequences {ap} and {bp} that is a
    refinement of the pairing {cdold->pv}.
    
    The new path is constrained to use only valid rungs and valid
    steps.
    
    The procedure fills the dynamic programming tableau{*tb}, which is
    expanded as needed. It also returns the last rung of the optimum
    pairing in {*goptp}, and its score in {*voptp}. The optimum
    pairing can be recovered by following the {prev} links in the
    tableau, starting with rung {*goptp}. 
    
    The parameters {delta,kapa,expand,shrink,maxUnp} are explained
    under {msm_cand_refine}. */

void msm_cand_refine_set_tableau_entry
  ( int r, /* {R}-coordinate of entry. */
    int s, /* {S}-coordinate of entry. */
    msm_seq_desc_t *ap,
    msm_seq_desc_t *bp,
    msm_pairing_t *prold,
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp,             /* Maximum unpaired datums between rungs. */
    msm_rung_step_score_proc_t *step_score, 
    bool_t may_begin_here,
    bool_t may_end_here,
    msm_rung_t *goptp,      /* IN/OUT: the current best choice for the last rung of pairing. */   
    double *voptp,          /* IN/OUT: score of that pairing, or {+INF}. */   
    msm_dyn_tableau_t *tb,  /* WORK: dynamic programming tableau. */
    int kp                  /* Hint for index of rung in {prold} whose {R}-coord is {~r}. */
  );
  /* Fills entry of {tb} corresponding to the rung {g} whose {R,S}
    coordinates are {r,s}. Assumes that {g} is a valid rung, given
    that the two sequences are {ap,bp} and the old psiring is {prold}.
    Unless said otherwise, the parameters have the same meanings as
    for {msm_cand_refine_fill_tableau}.
    
    If {may_begin_here} is TRUE, the procedure assumes that the optimum
    pairing is allowed to begin with the rung {g}.  
    
    If {may_end_here} is TRUE, the procedure assumes that the optimum 
    pairing may end with the rung {g}.  In that case, if the best pairing
    that ends at {g} has a better score than {*voptp}, the procedure 
    sets {*goptp} to {g} and {*voptp} to that score. */

void msm_cand_refine_compute_r_range
  ( msm_pairing_t *p, 
    int delta, 
    int kappa, 
    int *oldminrp, 
    int *oldmaxrp, 
    int *minrp, 
    int *maxrp
  );
  /* Computes the minimum and maximum R-coordinates {*oldminrp,*oldmaxrp}
    of the rungs in pairing {p}, and the minimum and maximum {R}-coordinates
    {*minrp,*maxrp} of any valid rung, given the parameters
    {delta,kappa}. */
  
void msm_refine_compute_s_range
  ( msm_pairing_t *p,
    int r,
    int *kp,
    int delta,
    int kappa,
    int *minsp,
    int *maxsp
  );
  /* Computes the minimum and maximum S-coordinates {*minsp,*maxsp} of
    any valid rung with R-coordinate {r}, given the old pairing {p}
    and the parameters {delta,kappa}. The variable {*kp} should
    contain a hint for the index of the rung of {p} whose
    {R}-coordinate closest to {r}; its value is updated by the
    procedure. */

/* IMPLEMENTATIONS */

msm_cand_t msm_cand_refine
  ( msm_cand_t *cdold, 
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp, 
    msm_rung_step_score_proc_t *step_score,
    msm_dyn_tableau_t *tb 
  )
  { 
    bool_t debug = FALSE; 
    
    if (debug)
      { fprintf(stderr, "  ----------------------------------------------------------\n");
        fprintf(stderr, "    refining candidate\n");
        msm_cand_write(stderr, "    old: ", cdold, 5, 10, 6, "\n");
      }
    
    /* We can't handle circular pairings yet: */
    bool_t circp = msm_pairing_is_circular(cdold->pr);
    demand((! circp), "circular pairings not supported yet");
    
    msm_rung_t gopt;
    double vopt;
    msm_cand_refine_fill_tableau
      ( cdold, delta, kappa, expand, shrink, maxUnp,
        step_score, &gopt, &vopt, tb
      );

    /* Get number of rungs {nr} in optimum pairing: */
    int nr = 0;
    msm_rung_t g = gopt;
    while (! msm_rung_is_none(g)) 
      { nr++; 
        msm_dyn_entry_t *ep = msm_dyn_tableau_get_entry_address(tb, g);
        g = (ep == NULL ? msm_rung_none : ep->prev);
      }
    /* Now collect the optimum pairing: */
    msm_rung_vec_t gvnew = msm_rung_vec_new(nr);
    g = gopt; int ip = nr;
    while (! msm_rung_is_none(g)) 
      { /* fprintf(stderr, "      ip = %3d g = (%5d,%5d)\n", ip, g.c[0], g.c[1]); */
        ip--; gvnew.e[ip] = g; 
        msm_dyn_entry_t *ep = msm_dyn_tableau_get_entry_address(tb, g);
        g = (ep == NULL ? msm_rung_none : ep->prev);
      }

    /* Create candidate: */
    msm_pairing_t *pr = msm_pairing_from_rung_vec(&gvnew, circp);
    msm_cand_t cdnew = msm_cand_from_pairing(&(cdold->seq[0]), &(cdold->seq[1]), pr, vopt); 

    if (debug)
      { fprintf(stderr, "\n");
        msm_cand_write(stderr, "    new: ", &cdnew, 5, 10, 6, "\n");
        fprintf(stderr, "  ----------------------------------------------------------\n");
      }
    return cdnew;
  }
  
void msm_cand_refine_fill_tableau
  ( msm_cand_t *cdold, 
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp,             /* Maximum unpaired datums between rungs. */
    msm_rung_step_score_proc_t *step_score, 
    msm_rung_t *goptp,      /* OUT: last rung of optimum pairing. */
    double *voptp,          /* OUT: score of optimum pairing. */
    msm_dyn_tableau_t *tb   /* WORK: dynamic programming tableau. */
  )
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "      filling matrices...\n"); }
    msm_pairing_t *prold = cdold->pr;
    
    demand(maxUnp >= 0, "bad {maxUnp}");
    
    /* We can't handle circular pairings yet: */
    demand((!msm_pairing_is_circular(prold)), "circular pairings not supported yet");
    
    /* Grab the two sequence descriptors: */
    msm_seq_desc_t *ap = &(cdold->seq[0]);
    msm_seq_desc_t *bp = &(cdold->seq[1]); 

    int oldminr, oldmaxr; /* Min and max R-coord of original rungs. */
    int minr, maxr; /* Min and max R-coord of any valid rung. */
    msm_cand_refine_compute_r_range(prold, delta, kappa, &oldminr, &oldmaxr, &minr, &maxr);

    /* Compute max differente {maxds} between two S-coordinates of valid rungs: */
    int maxds = 2*(kappa + 2*delta);
    /* Make sure that the tableau has all the required elements: */
    msm_dyn_tableau_resize(tb, minr, maxr, maxds);
    
    /* Fills the matrix within the band, in diagonal order.
      We process entries one diagonal at a time, in order of 
      increasing R-coordinate. */
    double vopt = -INF;   /* Score of the best path found so far. */
    msm_rung_t gopt = msm_rung_none;  /* Final rung {gopt} of that best path. */
    int r;
    int kp = 0; /* Approximate index in {p} of a rung with R-coord {r}; */
    for (r = minr; r <= maxr; r++)
      { /* Process tableau entries along the diagonal with R-coordinate {r}. */
        int mins, maxs; /* Range of S-coordinates for the R-coordinate {r}. */
        msm_refine_compute_s_range(prold, r, &kp, delta, kappa, &mins, &maxs);
        msm_dyn_tableau_set_s_range(tb, r, mins, maxs);

        /* Compute all elements in diagonal {r} of the tableau: */
        int s;
        for (s = mins; s <= maxs; s += 2)
          { assert((r + s) % 2 == 0);
            /* Check whether the path may begin or end here: */
            bool_t may_begin_here = (shrink || (r <= oldminr)) && (expand || (r >= oldminr));
            bool_t may_end_here = (shrink || (r >= oldmaxr)) && (expand || (r <= oldmaxr));
            /* Fill the tableau entry and update {gopt,vopt}: */
            msm_cand_refine_set_tableau_entry
              ( r, s,
                ap, bp, prold, 
                delta, kappa, expand, shrink, maxUnp, 
                step_score, 
                may_begin_here, may_end_here, 
                &gopt, &vopt, 
                tb, kp
              );
          }
      }
    
    /* Return last rung and score of optimum pairing: */
    if (debug) 
      { fprintf (stderr, "      vopt = %23.15f  gopt = (%5d,%5d)\n", vopt, gopt.c[0], gopt.c[1]); }
    (*goptp) = gopt; (*voptp) = vopt;
  }

void msm_cand_refine_set_tableau_entry
  ( int r, /* {R}-coordinate of entry. */
    int s, /* {S}-coordinate of entry. */
    msm_seq_desc_t *ap,
    msm_seq_desc_t *bp,
    msm_pairing_t *prold,
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp,             /* Maximum unpaired datums between rungs. */
    msm_rung_step_score_proc_t *step_score, 
    bool_t may_begin_here,
    bool_t may_end_here,
    msm_rung_t *goptp,         
    double *voptp,             
    msm_dyn_tableau_t *tb,  /* WORK: dynamic programming tableau. */
    int kp                  /* Hint for index of rung in {prold} whose {R}-coord is {~r}. */
  )
  {
    /* Compute the X and Y coordinates {ix,iy} of entry: */
    int ix = (r + s) / 2;
    int iy = (r - s) / 2;
    /* fprintf(stderr, "      ix = %5d iy = %5d\n", ix, iy); */
    msm_rung_t g = (msm_rung_t){{ix, iy}};
    msm_dyn_entry_t *eg = msm_dyn_tableau_get_entry_address(tb, g);
    assert(eg != NULL);
    bool_t ix_inside = ((ap->circ) || ((ix >= 0) && (ix < ap->npos)));
    bool_t iy_inside = ((bp->circ) || ((iy >= 0) && (iy < bp->npos)));
    if ((! ix_inside) || (! iy_inside))
      { /* Rung {g} is not valid, set tableau entry to null values: */
        eg->score = +INF; eg->prev = msm_rung_none;
      }
    else
      { /* Compute the maximum-score path ending at rung {g}. */

        double vsel = -INF;  /* Score of best path that ends at {g}. */
        msm_rung_t fsel;     /* Previous rung on that path. */

        auto void try_step(msm_rung_t f);
          /* Computes the score of the best pairing that ends
            with rung {f}, plus the cost of the step {f->g}.
            If the step is valid, updates {fsel,vsel}. */

        void try_step(msm_rung_t f)
          { 
            double vtof; /* Value of best path to {f}. */
            if (msm_rung_is_none(f))
              { vtof = 0; }
            else
              { msm_dyn_entry_t *ef = msm_dyn_tableau_get_entry_address(tb, f);
                vtof = (ef == NULL ? -INF : ef->score);
              }
            if (fabs(vtof) != INF)
              { /* Rung {f} is valid: */
                double vfg = step_score(ap, bp, f, g);
                /* fprintf(stderr, "\n"); */
                /* fprintf(stderr, "  step from (%d,%d)", f->c[0], f->c[1]); */
                /* fprintf(stderr, " to (%d,%d) ", g.c[0], g.c[1]); */
                /* fprintf(stderr, " score %9.6f + %9.6f = %9.6f", vtof, vfg, vf+vfg); */
                /* fprintf(stderr, "\n"); */
                double v = vtof + vfg;
                bool_t better;
                if (v > vsel)
                  { better = TRUE; }
                else if (v == vsel)
                  { better = msm_rung_step_break_tie(f, g, fsel, g); }
                else
                  { better = FALSE; }
                if (better) { vsel = v; fsel = f; }
              }
          }

        /* Choose between all valid maximum-score paths ending
          at rung {g}. Namely, consider coming from rung
          {msm_rung_none} (path starts at {g}) and from every
          previous valid rung {f} such that (1) the step
          {f->g} is strictly increasing, and (2) the step
          leaves at most {maxUnp} unpaired datums. */

        int minrf = r - (maxUnp + 2); /* Min {msm_rung_R(f)}. */
        int maxrf = r - 2;            /* Max {msm_rung_R(f)}. */
        int rf; /* R-coord of rung {f}. */
        int kpf = kp; /* Approximate index in {p} of a rung with R-coord {rf}; */
        for (rf = maxrf; rf >= minrf; rf--)
          { int minsf, maxsf; /* Range of S-coords for rung {f}. */
            /* Get range of S-coords in tableau for R-coord {rf}: */
            msm_refine_compute_s_range(prold, rf, &kpf, delta, kappa, &minsf, &maxsf);
            /* Scan S-coords of valid {f} rungs: */
            int sf;  /* S-coord of rung {f}. */
            for (sf = minsf; sf <= maxsf; sf += 2)
              { /* Compute rung {f} from {rf, sf}: */ 
                assert((rf + sf) % 2 == 0);
                int ixf = (rf + sf)/2;
                int iyf = (rf - sf)/2;
                if ((ixf < ix) && (iyf < iy))
                  { try_step((msm_rung_t){{ ixf, iyf }}); }
              }
          }
        if (may_begin_here) 
          { /* Consider starting here: */
            try_step(msm_rung_none);
          }

        /* Save the optimum path and the score in the matrix entry for {g}: */
        eg->score = vsel; eg->prev = fsel;
        /* fprintf(stderr, "      prev[%5d,%5d]", g.c[0], g.c[1]); */
        /* fprintf(stderr, " = (%5d,%5d)\n", gsel.c[0], gsel.c[1]); */

        if (may_end_here)
          { /* Update the best final rung {gopt} and score {vopt}: */
            if (vsel > (*voptp)) { (*voptp) = vsel; (*goptp) = g; }
          }
      }
  }

void msm_cand_refine_compute_r_range
  ( msm_pairing_t *p, 
    int delta, 
    int kappa, 
    int *oldminrp, 
    int *oldmaxrp, 
    int *minrp, 
    int *maxrp
  )
  {
    bool_t debug = FALSE;
    
    /* Get the first rung {glo} of pairing: */
    int nr = msm_pairing_fund_rung_count(p);
    msm_rung_t glo = msm_pairing_get_rung(p, 0);
    assert(! msm_pairing_is_circular(p));
    if (debug) { fprintf(stderr, "      lo rung = %5d %5d\n", glo.c[0], glo.c[1]); }
    
    /* Get the last rung {ghi} of pairing: */
    msm_rung_t ghi = msm_pairing_get_rung(p, nr-1);
    if (debug) { fprintf(stderr, "      hi rung = %5d %5d\n", ghi.c[0], ghi.c[1]); }
    
    /* Compute range {oldminr,oldmaxr} of R-coords of rungs in pairing {p}: */
    int oldminr = glo.c[0] + glo.c[1];
    int oldmaxr = ghi.c[0] + ghi.c[1];

    /* Compute minimum and maximum R-coord of any valid rung: */
    int minr = oldminr - 2*(delta + kappa);
    int maxr = oldmaxr + 2*(delta + kappa);

    /* Return results: */
    (*oldminrp) = oldminr;
    (*oldmaxrp) = oldmaxr;
    (*minrp) = minr;
    (*maxrp) = maxr;
  }
  
void msm_refine_compute_s_range
  ( msm_pairing_t *p,
    int r,
    int *kp,
    int delta,
    int kappa,
    int *minsp,
    int *maxsp
  )
  {
    /* Get num of rungs {nr} in given pairing: */
    int nr = msm_pairing_fund_rung_count(p);
    assert(! msm_pairing_is_circular(p));
    
    auto int rung_r(int k);
      /* Returns the R-coordinate ofrung {k} of pairing {p}. */
      
    auto int rung_s(int k);
      /* Returns the S-coordinate of rung {k} of pairing {p}. */
      
    int rung_r(int k)
      { msm_rung_t g = msm_pairing_get_rung(p, k);
        return msm_rung_R(g);
      }
    
    int rung_s(int k)
      { msm_rung_t g = msm_pairing_get_rung(p, k);
        return msm_rung_S(g);
      }

    /* Get the R- and S-coords {lor,los} of the first rung: */
    int lor = rung_r(0), los = rung_s(0);
    
    /* Get the R- and S-coords {hir,his} of the last rung */
    int hir = rung_r(nr-1), his = rung_s(nr-1);

    /* If {lor <= r <= hir}, find rung indices {kprev} and {knext},
      consecutive or equal, so that the step {p[kprev] -> p[knext]}
      crosses the sweepline; namely, such that {rung_r(kprev) <= r <=
      rung_r(knext)}. If {r < lor} we set {kprev} to {-1},
      if {r > hir} we set {knext} to {nr}; otherwise 
      we will have {0 <= kprev <= knext < nr}. */
    int kprev, knext;
    if (r < lor)
      { /* Sweepline lies before rung 0: */ 
        kprev = -1; knext = 0;
      }
    else if (r > hir)
      { /* Sweepline lies after rung {nr-1}: */
        kprev = nr-1; knext = nr;
      }
    else
      { /* Sweepline crosses the pairing {p}, must find {kprev,knext}. */
        /* Start the search for {kprev} with the given rung index {*kp}: */
        kprev = *kp; 
        if (kprev < 0) { kprev = 0; }
        if (kprev >= nr) { kprev = nr-1; }
        while ((kprev > 0) && (rung_r(kprev) > r)) { kprev--; }
        kprev++;
        while ((kprev < nr) && (rung_r(kprev) <= r)) { kprev++; }
        kprev--;
        assert((kprev >= 0) && (kprev < nr));
        
        /* Find {knext}: */
        knext = kprev;
        while ((knext < nr) && (rung_r(knext) < r)) { knext++; }
        assert(knext < nr);
      }
    /* Update approx index {*kp} for next call: */
    (*kp) = knext;

    /* Compute center {meds} and radius {rads} of the S-coord range for this R-coord: */
    int meds, rads;
    if (r <= lor - kappa)
      { meds = los; rads = 2*(kappa + delta) - (lor - r); }
    else if (r >= hir + kappa)
      { meds = his; rads = 2*(kappa + delta) - (r - hir); }
    else if (r < lor) 
      { meds = los; rads = lor - r + 2*delta; }
    else if (r > hir) 
      { meds = his; rads = r - hir + 2*delta; }
    else 
      { int sprev = rung_s(kprev), rprev = rung_r(kprev);
        int snext = rung_s(knext), rnext = rung_r(knext);
        assert((rprev <= r) && (r <= rnext));
        if (r == rprev) 
          { meds = sprev; rads = 2*delta; }
        else if (r == rnext) 
          { meds = snext; rads = 2*delta; }
        else 
          { /* Compute {meds} by linear interpolation: */
            int ds = snext - sprev;
            int dr = rnext - rprev;
            meds = sprev + ds*(r - rprev)/dr;
            /* The parity of {s} must agree with that of {r}: */
            if ((r + meds) % 2 == 0)
              { /* There is a {delta}-square centered here: */
                rads = 2*delta;
              }
            else
              { /* The nearest {delta}-square is centered at {r±1}: */
                rads = 2*delta - 1;
              }
          }
      }

    /* Compute range {mins..maxs} of the {s} coord, store in tableau: */
    int mins = meds - rads;
    int maxs = meds + rads;
    /* fprintf(stderr, "      r = %5d  mins = %5d maxs = %5d", r, mins, maxs); */
    /* fprintf(stderr, " (x,y)min = (%5d,%5d)", (r + mins)/2, (r - mins)/2); */
    /* fprintf(stderr, " (x,y)min = (%5d,%5d)", (r + maxs)/2, (r - maxs)/2); */
    /* fprintf(stderr, "\n"); */
    
    /* Return results: */
    (*minsp) = mins; (*maxsp) = maxs;
  }
