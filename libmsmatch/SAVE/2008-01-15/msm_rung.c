/* See msm_rung.h */
/* Last edited on 2008-01-12 08:27:04 by stolfi */

#define msm_rung_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <sign.h>
#include <fget.h>
#include <nget.h>
#include <jsrandom.h>

#include <msm_basic.h>
#include <msm_seq_desc.h>

#include <msm_rung.h>

sign_t msm_rung_break_tie( msm_rung_t *ga, msm_rung_t *gb);
  /* Returns 0 if both rungs are {msm_rung_none}, or both are different
    from {msm_rung_none}.  Returns {-1} if only {ga} is {msm_rung_none}.
    Returns {+1} if only {gb} is {msm_rung_none}. */

bool_t msm_rung_equal(msm_rung_t *ga, msm_rung_t *gb)
  { return ((ga->c[0] == gb->c[0]) && (ga->c[1] == gb->c[1])); }

bool_t msm_rung_is_none(msm_rung_t *g)
  { return ((g->c[0] == msm_rung_none.c[0]) && (g->c[1] == msm_rung_none.c[1])); }

msm_rung_t msm_rung_throw(msm_rung_t *glo, msm_rung_t *ghi, bool_t diagonal)
  { int dx = ghi->c[0] - glo->c[0];
    int dy = ghi->c[1] - glo->c[1];
    msm_rung_t g = *glo;
    if (diagonal)
      { /* Pick a point on the diagonal: */
        double r = drandom();
        g.c[0] += (int)floor(r*(dx + 1));
        g.c[1] += (int)floor(r*(dy + 1));
      }
    else
      { /* Pick a point anywhere in the rectangle: */
        g.c[0] += abrandom(0, dx); 
        g.c[1] += abrandom(0, dy);
      }
    return g;
  }

msm_rung_t msm_rung_throw_middle(msm_rung_t *glo, msm_rung_t *ghi, double skipProb)
  { /* Number of steps on each side: */
    int DX = ghi->c[0] - glo->c[0];
    int DY = ghi->c[1] - glo->c[1];
    demand(DX >= 2, "step X displacement is too small");
    demand(DY >= 2, "step Y displacement is too small");
    int nsmax = (DX < DY ? DX : DY); /* Max steps from {glo} to {ghi}. */
    int cs = (nsmax + abrandom(0,1))/2; /* Half the steps, rounded randomly */
    affirm((cs > 0) && (cs < nsmax), "bug");
    msm_rung_t g = *glo;
    if (DX > DY)
      { int dy = msm_choose(1, DY-1, skipProb);
        assert((dy > 0) && (dy < DY));
        g.c[1] += dy;
        g.c[0] += 2*cs - dy;
        if (g.c[0] <= glo->c[0]) { g.c[0]++; }
        if (g.c[0] >= ghi->c[0]) { g.c[0]--; }
      }
    else
      { int dx = msm_choose(1, DX-1, skipProb);
        assert((dx > 0) && (dx < DX));
        g.c[0] += dx;
        g.c[1] += 2*cs - dx; 
        if (g.c[1] <= glo->c[1]) { g.c[1]++; }
        if (g.c[1] >= ghi->c[1]) { g.c[1]--; }
      }
    /* Consistency check: */
    affirm(msm_rung_X(*glo) < msm_rung_X(g), "bug");
    affirm(msm_rung_X(g) < msm_rung_X(*ghi), "bug");
    affirm(msm_rung_Y(*glo) < msm_rung_Y(g), "bug");
    affirm(msm_rung_Y(g) < msm_rung_Y(*ghi), "bug");
    return g;
  }

msm_rung_vec_t msm_rung_vec_throw(msm_rung_t *gini, msm_rung_t *gfin, bool_t atomic, double skipProb)
  { /* Compute increments {DX,DY} in each coordinate: */
    int DX = gfin->c[0] - gini->c[0]; /* Increment in X. */
    int DY = gfin->c[1] - gini->c[1]; /* Increment in Y. */
    /* Allocate the rung vector: */
    msm_rung_vec_t gv = msm_rung_vec_new(DX + DY + 1);
    int ng = 0;
    /* Store initial rung: */
    gv.el[ng] = *gini; ng++;

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
        int dx = gf.c[0] - gi.c[0];
        int dy = gf.c[1] - gi.c[1];
        /* Check for empty case: */
        if ((dx == 0) && (dy == 0)) { /* Nothing to do. */ return; }
        /* Require strictly increasing: */
        assert((dx > 0) && (dy > 0)); 
        /* Compute the max number of steps {nsmax}: */
        int nsmax = (dx < dy ? dx : dy); /* Max steps from {gi} to {gf}. */
        assert(nsmax > 0); 
        if (nsmax == 1)
          { /* Just add the last rung {gf}: */
            msm_rung_vec_expand(&gv, ng);
            gv.el[ng] = gf; ng++;
          }
        else
          { /* Select a middle rung {gm} strictly between {gi} and {gm}: */
            msm_rung_t gm = msm_rung_throw_middle(&gi, &gf, skipProb);
            /* Add rungs from {gi} (excl) up to {gm} (incl): */
            append_rungs(gi, gm);
            /* Add rungs from {gm} (excl) up to {gf} (incl): */
            append_rungs(gm, gf);
          }
      }

    append_rungs(*gini, *gfin);
    msm_rung_vec_trim(&gv, ng);
    return gv;
  }

bool_t msm_rung_step_is_increasing(msm_rung_t *g0, msm_rung_t *g1, bool_t die)
  { int dx = g1->c[0] - g0->c[0];
    int dy = g1->c[1] - g0->c[1];
    if (dx < 1) { fail_test(die, "X step is too short"); }
    if (dy < 1) { fail_test(die, "Y step is too short"); }
    return TRUE;
  }

bool_t msm_rung_step_is_atomic(msm_rung_t *g0, msm_rung_t *g1, bool_t die)
  { int dx = g1->c[0] - g0->c[0];
    int dy = g1->c[1] - g0->c[1];
    if ((abs(dx) != 1) && (abs(dy) != 1)) { fail_test(die, "step is not atomic"); }
    return TRUE;
  }

bool_t msm_rung_step_is_perfect(msm_rung_t *g0, msm_rung_t *g1, bool_t die)
  { int dx = g1->c[0] - g0->c[0];
    int dy = g1->c[1] - g0->c[1];
    if ((dx != 1) || (dy != 1)) { fail_test(die, "step is not perfect"); }
    return TRUE;
  }

sign_t msm_rung_break_tie(msm_rung_t *ga, msm_rung_t *gb)
  { bool_t an = msm_rung_is_none(ga);
    bool_t bn = msm_rung_is_none(gb);
    if (an == bn)
      { return 0; }
    else if (an)
      { return -1; }
    else
      { return +1; }
  }

sign_t msm_rung_step_break_tie(msm_rung_t *a0, msm_rung_t *a1, msm_rung_t *b0, msm_rung_t *b1)
  { bool_t a0n = msm_rung_is_none(a0);
    bool_t a1n = msm_rung_is_none(a1);
    bool_t b0n = msm_rung_is_none(b0);
    bool_t b1n = msm_rung_is_none(b1);
        
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
        int dxa = a1->c[0] - a0->c[0];
        int dya = a1->c[1] - a0->c[1];
        int dxb = b1->c[0] - b0->c[0];
        int dyb = b1->c[1] - b0->c[1];
        /* Compare step lengths in L1 metric: */
        int la = abs(dxa) + abs(dya);
        int lb = abs(dxb) + abs(dyb);
        if (la < lb)
          { return +1; }
        else if (la > lb)
          { return -1; }
        /* Compare deviations from diagonal: */
        int va = abs(dxa - dya);
        int vb = abs(dxb - dyb);
        if (va < vb)
          { return +1; }
        else if (va > vb)
          { return -1; }
        /* Tie persists, give up: */
        return 0;
      }
  }

void msm_rung_interpolate(msm_rung_t g0, msm_rung_t g1, int *ngP, msm_rung_vec_t *gv)
  {
    int ng = (*ngP);
    /* Compute total displacement {sgnd[j]*absd[j]} in each axis {j}: */
    int absd[2], sgnd[2];
    int j;
    for (j = 0; j <= 1; j++)
      { int d = g1.c[j] - g0.c[j];
        absd[j] = abs(d);
        sgnd[j] = (d == 0 ? 0 : (d < 0 ? -1 : +1));
      }
    /* Find axes {jmin,jmax} with smallest and largest displacement: */
    int jmax = (absd[0] >= absd[1] ? 0 : 1);
    int jmin = 1 - jmax;
    demand(absd[jmin] >= 1, "step {g0-->g1} is stationary in X or Y");
    /* Step along axis {jmin}: */
    int dmin;
    for (dmin = 1; dmin <= absd[jmin]; dmin++)
      { /* Compute the relative position {dmax} along axis {jmax}: */
        int dmax = absd[jmax]*dmin/absd[jmin];
        /* Compute the absolute coordinates and pack into a rung {f}: */
        msm_rung_t f;
        f.c[jmin] = g0.c[jmin] + sgnd[jmin]*dmin;
        f.c[jmax] = g0.c[jmax] + sgnd[jmax]*dmax;
        /* Append the rung to {gv}: */
        msm_rung_vec_expand(gv, ng);
        gv->el[ng] = f; ng++;
      }
    (*ngP) = ng;
  }

msm_rung_vec_t msm_rung_vec_interpolate(msm_rung_vec_t *gv)
  {
    /* Get number of old rungs: */
    int ngold = gv->nel;
    /* Allocate the new rung vector {gvnew} (it may grow later): */
    msm_rung_vec_t gvnew = msm_rung_vec_new(ngold);
    if (ngold > 0)
      { /* The new rungs will be {gvnew[0..ngnew-1]}: */
        int ngnew = 0;
        /* Initialize the previous rung {g} with first rung of {gv}: */
        msm_rung_t g = gv->el[0];
        /* Insert the first rung in {gvnew}: */
        gvnew.el[0] = g; ngnew++;
        /* Process all steps of {gv}: */
        int i;
        for (i = 1; i < ngold; i++)
          { /* Get final rung {h} of this step: */
            msm_rung_t h = gv->el[i];
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

#define msm_rung_huge (msm_rung_t){{ INT_MAX, INT_MAX }}

msm_rung_vec_t msm_rung_vec_join(msm_rung_vec_t *gva, int ja, msm_rung_vec_t *gvb, int jb)
  { /* Allocate the result rung vector, with guessed size: */
    int ngmax = (gva->nel < gvb->nel ? gva->nel : gvb->nel);
    msm_rung_vec_t gv = msm_rung_vec_new(ngmax);
    /* Scan the two rung vectors and collect coincidences: */
    int ng = 0; /* Output rungs will be {gv.el[0..ng-1]}. */
    int ka = 0, kb = 0;  /* Indices into {gva,gvb}: */
    while ((ka < gva->nel) && (kb < gvb->nel))
      { msm_rung_t ga = (ka >= gva->nel ? msm_rung_huge : gva->el[ka]);
        msm_rung_t gb = (kb >= gvb->nel ? msm_rung_huge : gvb->el[kb]);
        if (ga.c[ja] < gb.c[jb]) 
          { ka++; }
        else if (ga.c[ja] > gb.c[jb])
          { kb++; }
        else
          { /* Found a matched pair of rungs: */
            msm_rung_t g = (msm_rung_t){{ ga.c[1-ja], gb.c[1-jb] }};
            if (ng > 0) { (void)msm_rung_step_is_increasing(&(gv.el[ng-1]), &g, /*die*/ TRUE); }
            msm_rung_vec_expand(&gv, ng);
            gv.el[ng] = g; ng++; 
            ka++; kb++; 
          }
      }
    msm_rung_vec_trim(&gv, ng);
    return gv;
  }

msm_rung_vec_t msm_rung_vec_map_coordinates(msm_rung_vec_t *gv, msm_rung_map_proc_t *map)
  { /* Get rung count {ng} and circularity {circp} of given rung_vec: */
    int ng = gv->nel;
    /* Allocate new rung vector: */
    msm_rung_vec_t gvnew = msm_rung_vec_new(ng);
    /* Map all rungs: */
    int i;
    for (i = 0; i < gvnew.nel; i++)
      { gvnew.el[i] = map(gv->el[i]); }
    /* Make rung_vec from rung vector: */
    return gvnew;
  }

sign_t msm_rung_lex_compare(msm_rung_t *ga, msm_rung_t *gb)
  { if (ga->c[0] < gb->c[0])
      { return NEG; }
    else if (ga->c[0] > gb->c[0])
      { return POS; }
    else if (ga->c[1] < gb->c[1])
      { return NEG; }
    else if (ga->c[1] > gb->c[1])
      { return POS; }
    else 
      { return ZER; }
  }

sign_t msm_rung_strict_compare(msm_rung_t *ga, msm_rung_t *gb)
  { if ((ga->c[0] == gb->c[0]) && (ga->c[1] == gb->c[1]))
      { return ZER; }
    else if ((ga->c[0] <= gb->c[0]) && (ga->c[1] <= gb->c[1]))
      { return NEG; }
    else if ((ga->c[0] >= gb->c[0]) && (ga->c[1] >= gb->c[1]))
      { return POS; }
    else 
      { return ZER; }
  }

#define msm_rung_MAXREP 3
 /* Maximum number of explicit character replications in {msm_rung_step_write}. */

void msm_rung_step_write(FILE *wr, msm_rung_t *g0, msm_rung_t *g1)
  {
    auto void wrstep(char cha, int rep, char chb);
      /* Prints a step as "{cha}{rep}{chb}" to {wr}.
      
        However, if {cha} is '\000', omits the "{cha}" part; if {rep}
        is 0, omits the "{cha}{rep}" part. Also, if {cha} is not
        '\000' and {rep} is in {1..msm_rung_MAXREP}, omits the "{rep}"
        part and repeats the "{cha}" part {rep} times. */
        
    void wrstep(char cha, int rep, char chb)
      { if (rep == 0) 
          { /* Omit {cha} and {rep}. */ }
        else if (cha == '\000') 
          { fprintf(wr, "%d", rep); }
        else if ((rep >= 1) && (rep <= msm_rung_MAXREP))
          { while (rep > 0) { fputc(cha, wr); rep--; } }
        else
          { fputc(cha, wr); fprintf(wr, "%d", rep); }
        fputc(chb, wr);
      }
    
    int dx = g1->c[0] - g0->c[0];
    int dy = g1->c[1] - g0->c[1];
    
    if ((dx == dy) && (dx > 0))
      { wrstep(':', dx-1, '|'); }
    else if (dx == 1)
      { wrstep('.', dy-1, '|'); }
    else if (dy == 1)
      { wrstep('\'', dx-1, '|'); }
    else if (dx == 0)
      { wrstep('\000', dy-1, '\\'); }
    else if (dy == 0)
      { wrstep('\000', dx-1, '/'); }
    else
      { /* Generic step notation: */
        fprintf(wr, "%d,%d|", dx, dy);
      }
  }

msm_rung_t msm_rung_step_read(FILE *rd, msm_rung_t g)
  { int ch = fgetc(rd);
    demand(ch != EOF, "unexpected EOF");
    char cha = '\000';     /* Leading char, or '\000' if none. */
    int rep = INT_MAX;  /* Repeat count of {cha}, or {INT_MAX} if none. */
    /* Read leading char {cha}, if any, and get its repeat count {rep}: */
    if ((ch == '\'') || (ch == '.') || (ch == ':'))
      { /* Gobble up the leading char {cha}, count repeats: */
        cha = ch; rep = 0;
        do 
          { rep++; ch = fgetc(rd); demand(ch != EOF, "unexpected EOF"); }
        while (ch == cha);
      }
    /* Read one integer into {rep}, or two integers into {dx,dy}: */
    int dx = INT_MAX, dy = INT_MAX; /* Explicit increment pair, or {INT_MAX} if none. */
    if (((ch >= '0') && (ch <= '9')) || (ch == '+') || (ch == '-'))
      { /* Read one integer into {rep}, or two integers into {dx,dy}: */
        demand((rep == INT_MAX) || (rep <= 1), "replicated {cha} with explicit {rep} or {dx,dy}"); 
        ungetc(ch, rd);
        int n;
        int ok = fscanf(rd, "%d", &n); 
        demand(ok != -1, "unexpected EOF");
        demand(ok == 1, "format error");
        ch = fgetc(rd); demand(ch != EOF, "unexpected EOF");
        if (ch == ',')
          { demand(cha == '\000', "explicit increments with non-null {cha}");
            dx = n; rep = INT_MAX;
            ok = fscanf(rd, "%d", &dy); 
            demand(ok != -1, "unexpected EOF");
            demand(ok == 1, "format error");
            ch = fgetc(rd); demand(ch != EOF, "unexpected EOF");
          }
        else
          { rep = n; }
      }
    /* At this point {ch} must be the main rung character ('\\', '/', or '|'): */
    if (ch == '/')
      { demand(cha == '\000', "invalid prefix char with \"/\" step");
        demand((dx == INT_MAX) && (dy == INT_MAX), "invalid \"/\" step");
        g.c[0] += (rep == INT_MAX ? 1 : rep);
      }
    else if (ch == '\\')
      { demand(cha == '\000', "invalid prefix char with \"\\\" step");
        demand((dx == INT_MAX) && (dy == INT_MAX), "invalid \"\\\" step");
        g.c[1] += (rep == INT_MAX ? 1 : rep);
      }
    else if (ch == '|')
      { if (cha == '\000')
          { demand(rep == INT_MAX, "invalid {rep} without {cha} and with \"|\""); 
            g.c[0] += (dx == INT_MAX ? 1 : dx);
            g.c[1] += (dy == INT_MAX ? 1 : dy);
          }
        else if (cha == ':')
          { demand((dx == INT_MAX) && (dy == INT_MAX), "invalid {dx,dy} with \":|\""); 
            g.c[0] += 1 + (rep == INT_MAX ? 1 : rep);
            g.c[1] += 1 + (rep == INT_MAX ? 1 : rep);
          }
        else if (cha == '\'')
          { demand((dx == INT_MAX) && (dy == INT_MAX), "invalid {dx,dy} with \"'|\""); 
            g.c[0] += 1 + (rep == INT_MAX ? 1 : rep);
            g.c[1] += 1;
          }
        else if (cha == '.')
          { demand((dx == INT_MAX) && (dy == INT_MAX), "invalid {dx,dy} with \".|\""); 
            g.c[0] += 1;
            g.c[1] += 1 + (rep == INT_MAX ? 1 : rep);
          }
        else
          { assert(FALSE); }
      }
    else 
      { demand(FALSE, "bad character"); }
    return g;
  }

msm_rung_vec_t msm_rung_vec_map_to_finer
  ( msm_rung_vec_t *gv, 
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
      
    return msm_rung_vec_map_coordinates(gv, &map_rung);
  }

msm_rung_vec_t msm_rung_vec_map_to_coarser
  ( msm_rung_vec_t *gv, 
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
      
    return msm_rung_vec_map_coordinates(gv, &map_rung);
  }
