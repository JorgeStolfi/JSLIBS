/* See msm_cand_refine.h */
/* Last edited on 2013-10-24 16:47:49 by stolfilocal */

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
#include <jsmath.h>

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
    int maxUnp,             /* Maximum unpaired datums between rungs. */
    msm_rung_step_score_proc_t *step_score, 
    msm_rung_t *goptp,      /* OUT: last rung of optimum pairing. */
    double *voptp,          /* OUT: score of optimum pairing. */
    msm_dyn_tableau_t *tb,  /* WORK: dynamic programming tableau. */
    int *n_steps,           /* Number of examined steps (in,out). */
    int *n_entries          /* Number of tableau entries that were computed (in,out). */
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
    bool_t debug,           /* TRUE to print diagnostics. */
    msm_rung_t *goptp,      /* IN/OUT: the current best choice for the last rung of pairing. */   
    double *voptp,          /* IN/OUT: score of that pairing, or {+INF}. */   
    msm_dyn_tableau_t *tb,  /* WORK: dynamic programming tableau. */
    int kp,                 /* Hint for index of rung in {prold} whose {R}-coord is {~r}. */
    int *n_steps,           /* Number of examined steps (in,out). */
    int *n_entries          /* Number of tableau entries that were computed. */
  );
  /* Fills entry of {tb} corresponding to the rung {g} whose {R,S}
    coordinates are {r,s}. Assumes that {g} is a valid rung for the tableau {tb},
    given that that the two sequences are {ap,bp} and the old pairing is {prold}.
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
    int n0, 
    int n1,
    int delta, 
    int kappa, 
    int *oldminrp, 
    int *oldmaxrp, 
    int *minrp, 
    int *maxrp
  );
  /* Computes the minimum and maximum R-coordinates {*oldminrp,*oldmaxrp}
    of the rungs in pairing {p}, and the minimum and maximum {R}-coordinates
    {*minrp,*maxrp} of any valid rung {(i0,i1}), given the parameters
    {delta,kappa}. Assumes that {i0} is restricted to the
    range {0..n0-1} and {i1} to the range {0..n1-1}. */
  
void msm_refine_compute_s_range
  ( msm_pairing_t *p,
    int n0, 
    int n1,
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
    contain a hint for the index of the rung {(i0,i1)} of {p} whose
    {R}-coordinate closest to {r}; its value is updated by the
    procedure. Assumes that {i0} is restricted to the
    range {0..n0-1} and {i1} to the range {0..n1-1}.
    May return an empty range. */

/* IMPLEMENTATIONS */

msm_cand_t msm_cand_refine
  ( msm_cand_t *cdold, 
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp, 
    msm_rung_step_score_proc_t *step_score,
    bool_t verbose,
    msm_dyn_tableau_t *tb,
    int *n_steps, 
    int *n_entries 
  )
  { 
    if (verbose)
      { fprintf(stderr, "  ----------------------------------------------------------\n");
        fprintf(stderr, "    refining candidate\n");
        msm_cand_write(stderr, "    old: ", cdold, 5, 10, 6, "\n", TRUE);
      }
    
    int n_steps_t = 0;
    int n_entries_t = 0;
    
    msm_rung_t gopt;
    double vopt;
    msm_cand_refine_fill_tableau
      ( cdold, delta, kappa, expand, shrink, maxUnp,
        step_score, &gopt, &vopt, tb,
        &n_steps_t, &n_entries_t
      );

    if(verbose)
      { fprintf
        ( stderr, "    %d entries, %.1f steps per entry\n", 
          n_entries_t, ((double)n_steps_t)/n_entries_t
        );
      }
    (*n_entries)+= n_entries_t;
    (*n_steps)+= n_steps_t;
    /* Get number of rungs {nr} in optimum pairing: */
    int nr = 0;
    msm_rung_t g = gopt;
    while (! msm_rung_is_none(&g)) 
      { nr++; 
        msm_dyn_entry_t *ep = msm_dyn_tableau_get_entry_address(tb, g);
        g = (ep == NULL ? msm_rung_none : ep->prev);
      }
    /* Now collect the optimum pairing: */
    msm_rung_vec_t gvnew = msm_rung_vec_new(nr);
    g = gopt; int ip = nr;
    while (! msm_rung_is_none(&g)) 
      { /* fprintf(stderr, "      ip = %3d g = (%5d,%5d)\n", ip, g.c[0], g.c[1]); */
        ip--; gvnew.e[ip] = g; 
        msm_dyn_entry_t *ep = msm_dyn_tableau_get_entry_address(tb, g);
        g = (ep == NULL ? msm_rung_none : ep->prev);
      }

    /* Create candidate: */
    msm_pairing_t *pr = msm_pairing_from_rung_vec(&gvnew);
    msm_cand_t cdnew = msm_cand_from_pairing(&(cdold->seq[0]), &(cdold->seq[1]), pr, vopt); 

    if (verbose)
      { fprintf(stderr, "\n");
        msm_cand_write(stderr, "    new: ", &cdnew, 5, 10, 6, "\n",TRUE);
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
    msm_dyn_tableau_t *tb,  /* WORK: dynamic programming tableau. */
    int *n_steps,           /* Number of examined steps (in,out). */
    int *n_entries          /* Number of tableau entries that were computed. */
  )
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "      filling matrices...\n"); }
    msm_pairing_t *prold = cdold->pr;
    
    demand(maxUnp >= 0, "bad {maxUnp}");
    
    /* Grab the two sequence descriptors: */
    msm_seq_desc_t *seq0 = &(cdold->seq[0]);
    msm_seq_desc_t *seq1 = &(cdold->seq[1]); 

    int n0 = seq0->size;
    int n1 = seq1->size;
    
    int oldminr, oldmaxr; /* Min and max R-coord of original rungs. */
    int minr, maxr; /* Min and max R-coord of any valid rung. */
    msm_cand_refine_compute_r_range(prold, n0, n1, delta, kappa, &oldminr, &oldmaxr, &minr, &maxr);

    /* Compute max differente {maxds} between two S-coordinates of valid rungs with same {R-coord}: */
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
        msm_refine_compute_s_range(prold, n0, n1, r, &kp, delta, kappa, &mins, &maxs);
        assert(mins <= maxs);
        msm_dyn_tableau_set_s_range(tb, r, mins, maxs);

        /* Compute all elements in diagonal {r} of the tableau: */
        int s;
        for (s = mins; s <= maxs; s += 2)
          { assert((r + s) % 2 == 0);
            /* Debugging decision: */
            int debug_drmin = 50;
            int debug_drmax = 100;
            bool_t debug_in_r_range = ((r - minr) >= debug_drmin) & ((r - minr) < debug_drmax);
            bool_t debug_in_s_range = (s == (mins + maxs)/2);
            bool_t debug_entry = debug_in_r_range & debug_in_s_range;
            /* Check whether the path may begin or end here: */
            bool_t may_begin_here = (shrink || (r <= oldminr)) && (expand || (r >= oldminr));
            bool_t may_end_here = (shrink || (r >= oldmaxr)) && (expand || (r <= oldmaxr));
            /* Fill the tableau entry and update {gopt,vopt}: */
            msm_cand_refine_set_tableau_entry
              ( r, s,
                seq0, seq1, prold, 
                delta, kappa, expand, shrink, maxUnp, 
                step_score, 
                may_begin_here, may_end_here, debug_entry,
                &gopt, &vopt, 
                tb, kp,
                n_steps, n_entries
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
    msm_seq_desc_t *seq0,
    msm_seq_desc_t *seq1,
    msm_pairing_t *prold,
    int delta,
    int kappa,
    bool_t expand,
    bool_t shrink,
    int maxUnp,             /* Maximum unpaired datums between rungs. */
    msm_rung_step_score_proc_t *step_score, 
    bool_t may_begin_here,
    bool_t may_end_here,
    bool_t debug,           /* TRUE to print diagnostics. */
    msm_rung_t *goptp,         
    double *voptp,             
    msm_dyn_tableau_t *tb,  /* WORK: dynamic programming tableau. */
    int kp,                 /* Hint for index of rung in {prold} whose {R}-coord is {~r}. */
    int *n_steps,           /* Number of examined steps (in,out). */
    int *n_entries          /* Number of tableau entries that were computed (in,out). */
  )
  {
    int n0 = seq0->size;
    int n1 = seq1->size;
    /* Compute the X and Y coordinates {i0,i1} of entry: */
    int i0 = (r + s) / 2;
    int i1 = (r - s) / 2;
    msm_rung_t g = (msm_rung_t){{ i0, i1 }};
    if (debug) { fprintf(stderr, "      filling entry [%d,%d]\n", i0, i1); }
    /* Try to get the tableau entry for {g}: */
    msm_dyn_entry_t *eg = msm_dyn_tableau_get_entry_address(tb, g);
    bool_t i0_inside = ((i0 >= 0) && (i0 < n0));
    bool_t i1_inside = ((i1 >= 0) && (i1 < n1));
    if ((eg == NULL) || (! i0_inside) || (! i1_inside))
      { /* Rung {g} is not valid: */
        if (eg != NULL) 
          { /* Set tableau entry to null values: */
            eg->score = -INF; eg->prev = msm_rung_none;
          }
        if (debug) { fprintf(stderr, "        entry is not valid\n"); }
      }
    else
      { /* Compute the maximum-score path ending at rung {g}. */

        double vsel = -INF;  /* Score of best path that ends at {g}. */
        msm_rung_t fsel = msm_rung_none;     /* Previous rung on that path. */

        auto void try_step(msm_rung_t f, int *kpfP);
          /* If {f} is inside the tableau and the step {f->g} is valid, 
            computes the score {vtofg} of the best pairing that ends
            with rung {f}, plus the cost of the step {f->g};
            if better than {vsel}, updates {fsel=f,vsel=vtofg}.
            The integer {*kpfP} is a hint to the position in {prold}
            of a rung with same R-coord as {f}. */

        void try_step(msm_rung_t f, int *kpfP)
          { double vtof; /* Value of best path to {f}. */
            if (msm_rung_is_none(&f))
              { if (debug) { fprintf(stderr, "        start at"); }
                vtof = 0;
              }
            else
              { if (debug) 
                  { fprintf(stderr, "        step from [%d,%d]", f.c[0], f.c[1]);
                    fprintf(stderr, " to [%d,%d]", g.c[0], g.c[1]);
                    fprintf(stderr, " (%d,%d)", g.c[0]-f.c[0], g.c[1]-f.c[1]);
                  }
                msm_dyn_entry_t *ef = msm_dyn_tableau_get_entry_address(tb, f);
                vtof = (ef == NULL ? -INF : ef->score);
              }
            if (vtof == -INF)
              { /* Rung {f} is outside tableau: */
                if (debug) { fprintf(stderr, " starts outside tableau"); }
              }
            else
              { /* Rung {f} is inside the tableau: */
                if (isnan(vtof) || (vtof == +INF))
                  { fprintf(stderr, "\n** f = (%d,%d) vtof = %9.4f\n", f.c[0], f.c[1], vtof);
                    assert(FALSE);
                  }
                double vfg = step_score(seq0, seq1, &f, &g);
                (*n_steps)++;
                double vtofg = vtof + vfg;
                if (debug) { fprintf(stderr, "  score %9.6f + %9.6f = %9.6f", vtof, vfg, vtofg); }
                bool_t better;
                if (vtofg > vsel)
                  { better = TRUE; }
                else if (vtofg == vsel)
                  { better = msm_rung_step_break_tie(f, g, fsel, g); }
                else
                  { better = FALSE; }
                if (better) { vsel = vtofg; fsel = f; }
              }
            if (debug) { fprintf(stderr, "\n"); }
          }
        
        /* Choose between all valid maximum-score paths ending
          at rung {g}. Namely, consider coming from rung
          {msm_rung_none} (path starts at {g}) and from every
          previous valid rung {f} such that (1) the step
          {f->g} is strictly increasing and atomic, and (2) the step
          leaves at most {maxUnp} unpaired datums all on 
          the same sequence. */
        
        int kpf = kp; /* Approximate index in {prold} of a rung with R-coord of {f}; */
        
        if (may_begin_here) 
          { /* Consider starting a pairing at {g}: */
            try_step(msm_rung_none, &kpf);
          }

        /* Consider a perfect step ending at {g}: */
        try_step((msm_rung_t){{ i0 - 1, i1 - 1 }}, &kpf);

        /* Consider the imperfect atomic steps ending at {g}: */
        int unp;
        for (unp = 1; unp <= maxUnp; unp++)
          { /* Consider the two {f} rungs that leave {unp} unpaired datums on each side: */
            int j;
            for (j = 0; j <= 1; j ++)
              { msm_rung_t f;
                f.c[j] = g.c[j] - 1 - unp;
                f.c[1-j] = g.c[1-j] - 1;
                try_step(f, &kpf);
              }
          }

        /* Save the optimum path and the score in the matrix entry for {g}: */
        eg->score = vsel; eg->prev = fsel;
        if (debug)
          { fprintf(stderr, "      best"); 
            if (msm_rung_is_none(&fsel))
              { fprintf(stderr, " start at [%d,%d]", g.c[0], g.c[1]); }
            else
              { fprintf(stderr, " from [%d,%d]", fsel.c[0], fsel.c[1]);
                fprintf(stderr, " (%d,%d)", g.c[0]-fsel.c[0], g.c[1]-fsel.c[1]);
              }
            fprintf(stderr, " score = %9.6f\n", vsel);
          }
        (*n_entries)++;
        
        if (may_end_here)
          { /* Update the overall best final rung {gopt} and best score {vopt}: */
            if (vsel > (*voptp)) { (*voptp) = vsel; (*goptp) = g; }
          }
      }
  }

void msm_cand_refine_compute_r_range
  ( msm_pairing_t *p, 
    int n0, 
    int n1,
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
    int nr = msm_pairing_num_rungs(p);
    msm_rung_t glo = msm_pairing_get_rung(p, 0);
    
    if (debug) { fprintf(stderr, "      lo rung = %5d %5d\n", glo.c[0], glo.c[1]); }
    
    /* Get the last rung {ghi} of pairing: */
    msm_rung_t ghi = msm_pairing_get_rung(p, nr-1);
    if (debug) { fprintf(stderr, "      hi rung = %5d %5d\n", ghi.c[0], ghi.c[1]); }
    
    /* Compute range {oldminr,oldmaxr} of R-coords of rungs in pairing {p}: */
    int oldminr = glo.c[0] + glo.c[1];
    int oldmaxr = ghi.c[0] + ghi.c[1];

    /* Compute minimum and maximum R-coord of any valid rung: */
    int min0 = (int)imax(glo.c[0] - kappa - delta,0);
    int min1 = (int)imax(glo.c[1] - kappa - delta,0);
    int max0 = (int)imin(ghi.c[0] + kappa + delta, n0 - 1);
    int max1 = (int)imin(ghi.c[1] + kappa + delta, n1 - 1);
    
    int minr = min0 + min1; 
    int maxr = max0 + max1;

      if(debug) 
        { fprintf(stderr,"range [%d..%d] expanded to [%d..to %d]\n", oldminr, oldmaxr, minr, maxr); }
    
    /* Return results: */
    (*oldminrp) = oldminr;
    (*oldmaxrp) = oldmaxr;
    (*minrp) = minr;
    (*maxrp) = maxr;
  }
  
void msm_refine_compute_s_range
  ( msm_pairing_t *p,
    int n0, 
    int n1,
    int r,
    int *kp,
    int delta,
    int kappa,
    int *minsp,
    int *maxsp
  )
  {
    /* Get num of rungs {nr} in given pairing: */
    int nr = msm_pairing_num_rungs(p);
    
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
    
    int min0 = (r + mins)/2;
    int max0 = (r + maxs)/2;
    int min1 = (r - maxs)/2;
    int max1 = (r - mins)/2;

    if(min0 < 0)   { int b = 0 - min0;       maxs = maxs - 2*b; min0 = min0 + b; max1 = max1 - b; }
    if(min1 < 0)   { int b = 0 - min1;       mins = mins + 2*b; max0 = max0 - b; min1 = min1 + b; }
    if(max0 >= n0) { int b = max0 - n0 + 1;  mins = mins + 2*b; max0 = max0 - b; min1 = min1 + b; }
    if(max1 >= n1) { int b = max1 - n1 + 1;  maxs = maxs - 2*b; min0 = min0 + b; max1 = max1 - b; }
    
    /* Return results: */
    (*minsp) = mins; (*maxsp) = maxs;
  }
