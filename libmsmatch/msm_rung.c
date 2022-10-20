/* See msm_rung.h */
/* Last edited on 2022-10-20 07:58:20 by stolfi */

#define msm_rung_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <sign.h>
#include <fget.h>
#include <nget.h>
#include <vec.h>
#include <jsrandom.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>

#include <msm_rung.h>

sign_t msm_rung_break_tie( msm_rung_t ga, msm_rung_t gb);
  /* Returns 0 if both rungs are {msm_rung_none}, or both are different
    from {msm_rung_none}.  Returns {-1} if only {ga} is {msm_rung_none}.
    Returns {+1} if only {gb} is {msm_rung_none}. */

bool_t msm_rung_equal(msm_rung_t ga, msm_rung_t gb)
  { return ((ga.c[0] == gb.c[0]) && (ga.c[1] == gb.c[1])); }

bool_t msm_rung_is_none(msm_rung_t *g)
  { return ((g->c[0] == msm_rung_none.c[0]) && (g->c[1] == msm_rung_none.c[1])); }

msm_rung_t msm_rung_throw(msm_rung_t glo, msm_rung_t ghi, bool_t diagonal)
  { int32_t d0 = ghi.c[0] - glo.c[0];
    int32_t d1 = ghi.c[1] - glo.c[1];
    msm_rung_t g = glo;
    if (diagonal)
      { /* Pick a point on the diagonal: */
        double r = drandom();
        g.c[0] += (int32_t)floor(r*(d0 + 1));
        g.c[1] += (int32_t)floor(r*(d1 + 1));
      }
    else
      { /* Pick a point anywhere in the rectangle: */
        g.c[0] += int32_abrandom(0, d0); 
        g.c[1] += int32_abrandom(0, d1);
      }
    return g;
  }

msm_rung_t msm_rung_throw_middle(msm_rung_t glo, msm_rung_t ghi, double skipProb)
  { /* Number of steps on each side: */
    int32_t D0 = ghi.c[0] - glo.c[0];
    int32_t D1 = ghi.c[1] - glo.c[1];
    demand(D0 >= 2, "step X displacement is too small");
    demand(D1 >= 2, "step Y displacement is too small");
    int32_t nsmax = (D0 < D1 ? D0 : D1); /* Max steps from {glo} to {ghi}. */
    int32_t cs = (nsmax + int32_abrandom(0,1))/2; /* Half the steps, rounded randomly */
    affirm((cs > 0) && (cs < nsmax), "bug");
    msm_rung_t g = glo;
    if (D0 > D1)
      { int32_t d1 = msm_choose(1, D1-1, skipProb);
        assert((d1 > 0) && (d1 < D1));
        g.c[1] += d1;
        g.c[0] += 2*cs - d1;
        if (g.c[0] <= glo.c[0]) { g.c[0]++; }
        if (g.c[0] >= ghi.c[0]) { g.c[0]--; }
      }
    else
      { int32_t d0 = msm_choose(1, D0-1, skipProb);
        assert((d0 > 0) && (d0 < D0));
        g.c[0] += d0;
        g.c[1] += 2*cs - d0; 
        if (g.c[1] <= glo.c[1]) { g.c[1]++; }
        if (g.c[1] >= ghi.c[1]) { g.c[1]--; }
      }
    /* Consistency check: */
    affirm(msm_rung_X(glo) < msm_rung_X(g), "bug");
    affirm(msm_rung_X(g) < msm_rung_X(ghi), "bug");
    affirm(msm_rung_Y(glo) < msm_rung_Y(g), "bug");
    affirm(msm_rung_Y(g) < msm_rung_Y(ghi), "bug");
    return g;
  }

msm_rung_vec_t msm_rung_vec_throw(msm_rung_t gini, msm_rung_t gfin, bool_t atomic, double skipProb)
  { /* Compute increments {D0,D1} in each coordinate: */
    int32_t D0 = gfin.c[0] - gini.c[0]; /* Increment in X. */
    int32_t D1 = gfin.c[1] - gini.c[1]; /* Increment in Y. */
    /* Allocate the rung vector: */
    msm_rung_vec_t gv = msm_rung_vec_new(D0 + D1 + 1);
    int32_t ng = 0;
    /* Store initial rung: */
    gv.e[ng] = gini; ng++;

    auto void append_rungs(msm_rung_t gi, msm_rung_t gf);
      /* Appends to {gv[0..ng-1]} a list of random rungs,
        interpolating between {gi} (exclusive) and {gf} (inclusive).
        All appended steps will be strictly increasing in both axes.
        If {atomic} is TRUE, the rungs will be atomic, too.
        
        The procedure assumes that the most recently appended rung
        ({gv[ng-1]}) is {gi}. If {gi==gf}, the procedure does nothing.
        Otherwise, if {gi-->gf} is atomic, the procedure merely appends
        {gf}. Otherwise, the step {gi-->gf} must be strictly increasing
        in both axes. */

    void append_rungs(msm_rung_t gi, msm_rung_t gf)
      { /* Steps on each side: */
        int32_t d0 = gf.c[0] - gi.c[0];
        int32_t d1 = gf.c[1] - gi.c[1];
        /* Check for empty case: */
        if ((d0 == 0) && (d1 == 0)) { /* Nothing to do. */ return; }
        /* Require strictly increasing: */
        assert((d0 > 0) && (d1 > 0)); 
        /* Compute the max number of steps {nsmax}: */
        int32_t nsmax = (d0 < d1 ? d0 : d1); /* Max steps from {gi} to {gf}. */
        assert(nsmax > 0); 
        if (nsmax == 1)
          { /* Just add the last rung {gf}: */
            msm_rung_vec_expand(&gv, ng);
            gv.e[ng] = gf; ng++;
          }
        else
          { /* Select a middle rung {gm} strictly between {gi} and {gm}: */
            msm_rung_t gm = msm_rung_throw_middle(gi, gf, skipProb);
            /* Add rungs from {gi} (excl) up to {gm} (incl): */
            append_rungs(gi, gm);
            /* Add rungs from {gm} (excl) up to {gf} (incl): */
            append_rungs(gm, gf);
          }
      }

    append_rungs(gini, gfin);
    msm_rung_vec_trim(&gv, ng);
    return gv;
  }

bool_t msm_rung_step_is_increasing(msm_rung_t g0, msm_rung_t g1, int32_t minIncrEach, int32_t minIncrSum, bool_t die)
  { int32_t d0 = g1.c[0] - g0.c[0];
    int32_t d1 = g1.c[1] - g0.c[1];
    if (d0 < minIncrEach) { fail_test(die, "X increment is too short"); }
    if (d1 < minIncrEach) { fail_test(die, "Y increment is too short"); }
    if (d0 + d1 < minIncrSum) { fail_test(die, "total X+Y increment are too short"); }
    return TRUE;
  }

bool_t msm_rung_step_is_atomic(msm_rung_t g0, msm_rung_t g1, bool_t die)
  { int32_t d0 = g1.c[0] - g0.c[0];
    int32_t d1 = g1.c[1] - g0.c[1];
    if ((abs(d0) != 1) && (abs(d1) != 1)) { fail_test(die, "step is not atomic"); }
    return TRUE;
  }

bool_t msm_rung_step_is_perfect(msm_rung_t g0, msm_rung_t g1, bool_t die)
  { int32_t d0 = g1.c[0] - g0.c[0];
    int32_t d1 = g1.c[1] - g0.c[1];
    if ((d0 != 1) || (d1 != 1)) { fail_test(die, "step is not perfect"); }
    return TRUE;
  }
  
int32_t msm_rung_step_span_increment(msm_rung_t g0, msm_rung_t g1, int32_t j)
  { if ((j == 0) || (j == 1))
      { return g1.c[j] - g0.c[j]; }
    else if (j == -1)
      { return (g1.c[0] - g0.c[0]) + (g1.c[1] - g0.c[1]); }
    else
      { demand(FALSE, "bad side"); return 0; }
  }

sign_t msm_rung_break_tie(msm_rung_t ga, msm_rung_t gb)
  { bool_t an = msm_rung_is_none(&ga);
    bool_t bn = msm_rung_is_none(&gb);
    if (an == bn)
      { return 0; }
    else if (an)
      { return -1; }
    else
      { return +1; }
  }

sign_t msm_rung_step_break_tie(msm_rung_t a0, msm_rung_t a1, msm_rung_t b0, msm_rung_t b1)
  { bool_t a0n = msm_rung_is_none(&a0);
    bool_t a1n = msm_rung_is_none(&a1);
    bool_t b0n = msm_rung_is_none(&b0);
    bool_t b1n = msm_rung_is_none(&b1);
        
    /* Check for rungs to/from {msm_rung_none}: */
    if ((a0n || a1n) && (b0n || b1n))
      { /* Both steps have a 'none' end-rung --- they are incomparable: */
        return 0; 
      }
    else if (a0n || a1n)
      { /* Step {a} (only) has a 'none' end-rung: */
        return -1;
      }
    else if (b0n || b1n)
      { /* Step {b} (only) has a 'none' end-rung: */
        return +1;
      }
    else
      { /* All end-rungs are valid. Compute the displacements: */
        int32_t d0a = a1.c[0] - a0.c[0];
        int32_t d1a = a1.c[1] - a0.c[1];
        int32_t d0b = b1.c[0] - b0.c[0];
        int32_t d1b = b1.c[1] - b0.c[1];
        /* Compare step lengths in L1 metric: */
        int32_t la = abs(d0a) + abs(d1a);
        int32_t lb = abs(d0b) + abs(d1b);
        if (la < lb)
          { return +1; }
        else if (la > lb)
          { return -1; }
        /* Compare deviations from diagonal: */
        int32_t va = abs(d0a - d1a);
        int32_t vb = abs(d0b - d1b);
        if (va < vb)
          { return +1; }
        else if (va > vb)
          { return -1; }
        /* Tie persists, give up: */
        return 0;
      }
  }

void msm_rung_interpolate(msm_rung_t g0, msm_rung_t g1, int32_t *ngP, msm_rung_vec_t *gv)
  {
    int32_t ng = (*ngP);
    /* Compute total displacement {sgnd[j]*absd[j]} in each axis {j}: */
    int32_t absd[2], sgnd[2];
    int32_t j;
    for (j = 0; j <= 1; j++)
      { int32_t d = g1.c[j] - g0.c[j];
        absd[j] = abs(d);
        sgnd[j] = (d == 0 ? 0 : (d < 0 ? -1 : +1));
      }
    /* Find axes {jmin,jmax} with smallest and largest displacement: */
    int32_t jmax = (absd[0] >= absd[1] ? 0 : 1);
    int32_t jmin = 1 - jmax;
    demand(absd[jmin] >= 1, "step {g0-->g1} is stationary in X or Y");
    /* Step along axis {jmin}: */
    int32_t dmin;
    for (dmin = 1; dmin <= absd[jmin]; dmin++)
      { /* Compute the relative position {dmax} along axis {jmax}: */
        int32_t dmax = absd[jmax]*dmin/absd[jmin];
        /* Compute the absolute coordinates and pack into a rung {f}: */
        msm_rung_t f;
        f.c[jmin] = g0.c[jmin] + sgnd[jmin]*dmin;
        f.c[jmax] = g0.c[jmax] + sgnd[jmax]*dmax;
        /* Append the rung to {gv}: */
        msm_rung_vec_expand(gv, ng);
        gv->e[ng] = f; ng++;
      }
    (*ngP) = ng;
  }

msm_rung_vec_t msm_rung_vec_interpolate(msm_rung_vec_t *gv)
  {
    /* Get number of old rungs: */
    int32_t ngold = gv->ne;
    /* Allocate the new rung vector {gvnew} (it may grow later): */
    msm_rung_vec_t gvnew = msm_rung_vec_new(ngold);
    if (ngold > 0)
      { /* The new rungs will be {gvnew[0..ngnew-1]}: */
        int32_t ngnew = 0;
        /* Initialize the previous rung {g} with first rung of {gv}: */
        msm_rung_t g = gv->e[0];
        /* Insert the first rung in {gvnew}: */
        gvnew.e[0] = g; ngnew++;
        /* Process all steps of {gv}: */
        int32_t i;
        for (i = 1; i < ngold; i++)
          { /* Get final rung {h} of this step: */
            msm_rung_t h = gv->e[i];
            /* Interpolate step {g->h} with atomic steps: */
            msm_rung_interpolate(g, h, &ngnew, &gvnew);
            /* Prepare for next step: */
            g = h;
          }
        /* Cut {gvnew} to the exact size: */
        msm_rung_vec_trim(&gvnew, ngnew);
      }
    return gvnew;
  }

msm_rung_vec_t msm_rung_vec_make_increasing(msm_rung_vec_t *gv, int32_t minIncrEach, int32_t minIncrSum)
  { int32_t ng = gv->ne;
    assert(ng > 0);
    
    /* Create the rung vector {gv_new} of result: */
    msm_rung_vec_t gv_new = msm_rung_vec_new(ng);
    int32_t ng_new = 0; /* The gathered rungs are {gv_new.e[0..ng_new-1]}. */
            
    /* Get first rung and last rungs {gini,gfin} of {gv}: */
    msm_rung_t gini = gv->e[0];
    msm_rung_t gfin = gv->e[ng-1];
    
    int32_t tinc0 = gfin.c[0] - gini.c[0]; /* Total increment on side 0. */
    int32_t tinc1 = gfin.c[1] - gini.c[1]; /* Total increment on side 1. */
    if ((tinc0 < minIncrEach) || (tinc1 < minIncrEach) || (tinc0+tinc1 < minIncrSum))
      { /* Result has a single rung: */
        msm_rung_t gmid = gv->e[ng/2];
        gv_new.e[0] = gmid;
        ng_new = 1;
      }
    else
      { /* Result has at least two rungs, {gini} and {gfin}. */
        /* Take {gini}: */
        gv_new.e[0] = gini;
        ng_new = 1; /* The gathered rungs are {gv_new.e[0..ng_new-1]}. */
        /* Gather the intermediate rungs that fit: */
        msm_rung_t gpre = gini; /* Last rung that was taken. */
        int32_t k = 1;
        while (k < ng - 1)
          { /* Get rung {k} of {gv}: */
            msm_rung_t gk = gv->e[k];
            bool_t take = TRUE; /* Should we take this rung? */
            /* Get increments {dpre[0..1],dpos[0..1]} between {gpre,gk,gfin}: */
            int32_t dpre[2], dpos[2];
            int32_t j;
            for (j = 0; j < 2; j++)
              { dpre[j] = gk.c[j] - gpre.c[j];
                dpos[j] = gfin.c[j] - gk.c[j];
                if( dpre[j] < 0 ){
                  fprintf(stderr,"Step %d - [%d,%d] -> [%d,%d]\n",k,gpre.c[0],gpre.c[1],gk.c[0],gk.c[1]);
                }
                demand(dpre[j] >= 0, "pairing is decreasing");
                /* Check whether the step lengths before and after on each seq are acceptable: */
                if ((dpre[j] < minIncrEach) || (dpos[j] < minIncrEach)) { take = FALSE; }
              }
            /* Check whether the total step lengths before and after on both seqs are acceptable: */
            if (dpre[0] + dpre[1] < minIncrSum) { take = FALSE; }
            if (dpos[0] + dpos[1] < minIncrSum) { take = FALSE; }
            /* So, shall we take it? */
            if (take)
              { msm_rung_vec_expand(&gv_new, ng_new); 
                gv_new.e[ng_new] = gk; ng_new++;
                gpre = gk;
              }
            /* In any case, go to next rung: */
            k++; 
          }
        /* Add the last rung: */
        msm_rung_vec_expand(&gv_new, ng_new); 
        gv_new.e[ng_new] = gfin; ng_new++;
      }
    /* Trim and return: */
    msm_rung_vec_trim(&gv_new, ng_new);
    return gv_new;
  }

#define msm_rung_huge (msm_rung_t){{ INT32_MAX, INT32_MAX }}

msm_rung_vec_t msm_rung_vec_join(msm_rung_vec_t *gva, int32_t ja, msm_rung_vec_t *gvb, int32_t jb)
  { /* Allocate the result rung vector, with guessed size: */
    int32_t ngmax = (gva->ne < gvb->ne ? gva->ne : gvb->ne);
    msm_rung_vec_t gv = msm_rung_vec_new(ngmax);
    /* Scan the two rung vectors and collect coincidences: */
    int32_t ng = 0; /* Output rungs will be {gv.e[0..ng-1]}. */
    int32_t ka = 0, kb = 0;  /* Indices into {gva,gvb}: */
    while ((ka < gva->ne) && (kb < gvb->ne))
      { msm_rung_t ga = (ka >= gva->ne ? msm_rung_huge : gva->e[ka]);
        msm_rung_t gb = (kb >= gvb->ne ? msm_rung_huge : gvb->e[kb]);
        if (ga.c[ja] < gb.c[jb]) 
          { ka++; }
        else if (ga.c[ja] > gb.c[jb])
          { kb++; }
        else
          { /* Found a matched pair of rungs: */
            msm_rung_t g = (msm_rung_t){{ ga.c[1-ja], gb.c[1-jb] }};
            if (ng > 0) { (void)msm_rung_step_is_increasing(gv.e[ng-1], g, 1, 1, /*die*/ TRUE); }
            msm_rung_vec_expand(&gv, ng);
            gv.e[ng] = g; ng++; 
            ka++; kb++; 
          }
      }
    msm_rung_vec_trim(&gv, ng);
    return gv;
  }

sign_t msm_rung_lex_compare(msm_rung_t ga, msm_rung_t gb)
  { if (ga.c[0] < gb.c[0])
      { return NEG; }
    else if (ga.c[0] > gb.c[0])
      { return POS; }
    else if (ga.c[1] < gb.c[1])
      { return NEG; }
    else if (ga.c[1] > gb.c[1])
      { return POS; }
    else 
      { return ZER; }
  }

sign_t msm_rung_strict_compare(msm_rung_t ga, msm_rung_t gb)
  { if ((ga.c[0] == gb.c[0]) && (ga.c[1] == gb.c[1]))
      { return ZER; }
    else if ((ga.c[0] <= gb.c[0]) && (ga.c[1] <= gb.c[1]))
      { return NEG; }
    else if ((ga.c[0] >= gb.c[0]) && (ga.c[1] >= gb.c[1]))
      { return POS; }
    else 
      { return ZER; }
  }

msm_rung_vec_t msm_rung_vec_map_gen(msm_rung_vec_t *gv, msm_rung_map_proc_t *map)
  { /* Get rung count {ng} of given rung_vec: */
    int32_t ng = gv->ne;
    /* Allocate new rung vector: */
    msm_rung_vec_t gvnew = msm_rung_vec_new(ng);
    /* Map all rungs: */
    int32_t i;
    for (i = 0; i < gvnew.ne; i++) { gvnew.e[i] = map(gv->e[i]); }
    /* Make rung_vec from rung vector: */
    return gvnew;
  }

msm_rung_vec_t msm_rung_vec_map
  ( msm_rung_vec_t *gv, 
    msm_seq_desc_t *s0_old,
    msm_seq_desc_t *s1_old,
    msm_seq_desc_t *s0_new,
    msm_seq_desc_t *s1_new
  )
  { 
    auto msm_rung_t map_rung(msm_rung_t g);
    
    msm_rung_t map_rung(msm_rung_t g)
      { double p0 = msm_seq_desc_map_index(g.c[0], s0_old, s0_new);
        double p1 = msm_seq_desc_map_index(g.c[1], s1_old, s1_new);
        int32_t i0 = (int32_t)msm_round(p0);
        int32_t i1 = (int32_t)msm_round(p1);
        return (msm_rung_t){{ i0, i1 }};
      }
      
    return msm_rung_vec_map_gen(gv, &map_rung);
  }

#define msm_rung_MAXD 3
  /* Maximum number of explicit ":", ".", or "'" 
    character replications in {msm_rung_step_write}. */

#define msm_rung_MAXREP 6
  /* Maximum number of explicit "|" character replications
    in {msm_rung_step_write}. */

void msm_rung_step_write(FILE *wr, int32_t d0, int32_t d1, int32_t rep)
  {
    auto void wrstep(char cha, int32_t d, char chb);
      /* Prints a step as "{cha}{d}{chb}" to {wr}.  Usually
        {cha} is ':', '.', or '\''; and {chb} 
        is  |', '/', or '\\'.
      
        However, if if {d} is 0, omits the "{cha}{d}" part;
        if {d} is in {1..msm_rung_MAXD}, omits the "{d}"
        part and repeats the "{cha}" part {d} times. */
        
    void wrstep(char cha, int32_t d, char chb)
      { if (d == 0) 
          { /* Omit {cha} and {d}. */ }
        else if ((d >= 1) && (d <= msm_rung_MAXD))
          { while (d > 0) { fputc(cha, wr); d--; } }
        else
          { fputc(cha, wr); fprintf(wr, "%d", d); }
        fputc(chb, wr);
      }
    
    demand((d0 >= 0) && (d1 >= 0), "step increments cannot be negative");
    demand(rep > 0, "step count must be positive");
    if ((d0 == d1) && (d0 > 0))
      { wrstep(':', d0-1, '|'); }
    else if (d0 == 1)
      { wrstep('.', d1-1, '|'); }
    else if (d1 == 1)
      { wrstep('\'', d0-1, '|'); }
    else if ((d0 == 0) && (d1 > 0))
      { wrstep('.', d1-1, '\\'); }
    else if ((d0 > 0) && (d1 == 0))
      { wrstep('\'', d0-1, '/'); }
    else
      { /* Generic step notation: */
        fprintf(wr, "(%d,%d)|", d0, d1);
      }
    if (rep > msm_rung_MAXREP)
      { /* Write num rungs in perfect segmetn; */
        fprintf(wr, "%d", rep); 
      }
    else if (rep > 1)
      { /* Print the extra rungs in perfect segment; */
        int32_t k;
        for (k = 1; k < rep; k++) { fprintf(wr, "|"); }
      }
  }

void msm_rung_step_read(FILE *rd, int32_t *d0P, int32_t *d1P, int32_t *repP)
  { int32_t ch = fgetc(rd); /* Next char in {rd}. */
    demand(ch != EOF, "unexpected EOF");
    /* Parse the increments {d0,d1}: */
    /* Interpret {d0} and {d1} as if next char was '|', fix later if not: */
    int32_t d0, d1;
    if (ch == '(')
      { /* Parse "{d0},{d1})|" (no '\\' or '/'): */
        ungetc(ch, rd);
        d0 = fget_int32(rd);
        fget_match(rd, ",");
        d1 = fget_int32(rd);
        fget_match(rd, ")");
        /* Get hold again of the next char: */
        ch = fgetc(rd); 
        demand(ch != EOF, "unexpected EOF");
      }
    else if ((ch == '\'') || (ch == '.') || (ch == ':'))
      { /* Save the leading char {cha}: */
        char cha = (char)ch;    /* Leading char. */
        /* Get the number of samples skipped {d} by counting instances of {cha}: */
        int32_t d = 0;
        do 
          { d++; 
            ch = fgetc(rd);
            demand(ch != EOF, "unexpected EOF");
          }
        while (ch == cha);
        /* Check for explicit sample skip count: */
        if ((d == 1) && ((ch >= '0') && (ch <= '9')))
          { /* There is an explicit sample skip count, parse it: */
            ungetc(ch, rd);
            d = fget_int32(rd);
            /* Get hold again of the next char: */
            ch = fgetc(rd); 
            demand(ch != EOF, "unexpected EOF");
          }
        if (cha == ':')
          { d0 = d + 1; d1 = d + 1; }
        else if (cha == '\'')
          { d0 = d + 1; d1 = 1; }
        else if (cha == '.')
          { d0 = 1; d1 = d + 1;}
        else
          { assert(FALSE); }
      }
    else
      { /* Neither "(d0,d1)" prefix nor [.':] prefixes: */
        d0 = 1; d1 = 1;
      }
    /* Now {ch} is the next char in {rd} and is not {EOF}. */
    
    /* Parse the rung char, fix {d0,d1}: */
    if (ch == '|')
      { /* Leave {d0,d1} alone. */ }
    else if (ch == '/')
      { demand(d1 == 1, "invalid use of \".\" with \"/\"");
        d1 = 0;
      }
    else if (ch == '\\')
      { demand(d0 == 1, "invalid use of \"'\" with \"\\\"");
        d0 = 0;
      }
    else
      { demand(FALSE, "invalid rung character"); } /*  */
    /* Now {ch} must have been consumed. */
    
    /* Parse the (numeric or explicit) count {rep} of rungs in the perfect segment: */
    int32_t rep = 1;
    ch = fgetc(rd);
    if (ch != EOF)
      { if ((ch >= '0') && (ch <= '9'))
          { /* Numeric rep field: */
            ungetc(ch, rd);
            rep = fget_int32(rd); 
            demand(rep >= 1, "field {rep} must be positive"); 
          }
        else
          { /* Count explicit '|' rung marks (no spaces allowed): */
            while (ch == '|')
              { rep++;
                ch = fgetc(rd);
              }
            if (ch != EOF) { ungetc(ch, rd); }
          }
      }
    
    /* Return results: */
    demand((d0 >= 0) && (d1 >= 0), "pairing is not monotonic");
    (*d0P) = d0;
    (*d1P) = d1;
    (*repP) = rep;
  }

vec_typeimpl(msm_rung_vec_t,msm_rung_vec,msm_rung_t);
