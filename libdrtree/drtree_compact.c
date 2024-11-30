/* See {drtree_compact.h} */
/* Last edited on 2023-06-24 11:08:17 by stolfi */

#define drtree_compact_C_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <cmp.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <in.h>
#include <vec.h> 

#include <drtree.h>
#include <drtree_compact.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

void drtree_compact_arrange
  ( int32_t ni, 
    drtree_node_t dt[], 
    int32_t tMin, 
    int32_t tMax, 
    int32_t rdr[],
    int32_t *ncols_P,
    int32_t *nrows_P
  )
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "    > %s\n", __FUNCTION__); }

    demand((ni >= 0) && (ni < drtree_indivs_MAX), "invalid {ni}");
    demand(tMin < tMax, "invalid {tMin..tMax}");

    /* Count children: */
    int32_t *nch = drtree_count_children(ni, dt);
    
    /* Allocate the cell occupancy array {occ}: */
    /* Use table {occ[row*ncols + col} to tell whether a cell is occupied:  */ 
    int32_t ncols = tMax - tMin + 1;
    int32_t nrows = ni; /* Tentative -- to be reduced later. */
    int32_t ncells = nrows*ncols; /* Number of cells in main plot area. */
    bool_t *occ = (bool_t *)notnull(malloc(ncells*sizeof(bool_t)), "no mem"); /* Cell occupancy */

    /* Table of displacements used by {find_free_row}: */
    int32_t nd = 3*nrows; /* Num of row displacements in table. */
    int32_t drow[nd];
    for (uint32_t kc = 0;  kc < nrows; kc++) 
     { /* Two attempts below, then one above: */
       assert(3*kc+2 < nd);
       drow[3*kc+0] = -2*kc-1; 
       drow[3*kc+1] = -2*kc-2; 
       drow[3*kc+2] = kc+1;
     } 
     assert(drow[3*nrows-1] >= +(nrows-1));
     assert(drow[3*nrows-2] <= -(nrows-1));

    /* Used by {assign_row}, see below. */
    double rChSum[ni];

    auto int32_t assign_row(int32_t iq);
      /* Assigns a row for the trace of an individual {iq}, so that the
        cells in columns corresponding to its life span are free,
        with one free cell to spare on each side.
        
        Must be called for {iq} only after it has been called for all
        nodes after node {iq}.
        
        The procedure assumes that, for all {jq} in {iq+1..ni-1}, {rdr[jq]} is
        is {-1} if individual {jq} is not visible, oherwise it is the 
        row that was assigned to {jq} in this iteration
        
        If the node {iq} is not visible, {rdr[iq]} must be {-1}, and the
        procedure returns {-1}. 
        
        If the node {iq} is visible, it assumes that {rdr[iq]} is the row
        that was assigned to {iq} in a previous iteration, or some
        arbitary row index in {0..nrows-1} for the first iteration.
        
        Moreover, for any non-null node {ip} in {0..ni-1} also
        assumes that {rChSum[ip]} is the sum of {rdr[ir]} for all
        children {ir} of {ip} that have already been assigned. In
        particular, when {ip = iq}, {rChSum} accounts for all
        {nch[iq]} children of {iq}.
        
        Then, if {iq} has any children, the ideal row for {iq} is taken
        to be {(rChSum[iq] + rdr[iq])/(nch[iq] + 1)}, rounded. If {iq}
        has no children, but has a parent {ip}, the ideal row is set to
        be {rdr[ip]}. If {iq} has no chidlren and no parent, its ideal
        row is row 0.
        
        In any case, the procedure will to assign {iq} to the row that
        has the required cells still free and is closest to that ideal
        row.  Then update {rChSum} of {iq}'s parent, if any. */

    auto void unassign_row(int32_t iq);
      /* If {rdr[iq]} is {-1} does nothing. Otherwise clears the cells
        currently used by {iq} and sets {rdr[iq]} to {-1}. */

    /* Assign a row {rdr[iq]} to every individual {iq}: */
    int32_t niters = 10; /* Iterations of the placement loop. */
    /* Pretend the initial placement was row 0: */
    for (uint32_t iq = 0;  iq < ni; iq++) { rdr[iq] = 0; }
    for (uint32_t it = 0;  it < niters; it++)
      { if (debug) { fprintf(stderr, "      --- iteration %d ---\n", it); }
        /* Assign rows in reverse chrono order (children before parent): */
        /* Clear {occ,rChSum}: */
        for (uint32_t k = 0;  k < ncells; k++) { occ[k] = FALSE; }
        for (uint32_t iq = 0;  iq < ni; iq++) { rChSum[iq] = 0.0; }
        for (int32_t iq = ni-1; iq >= 0; iq--) 
          { int32_t row = assign_row(iq);
            if (debug) { fprintf(stderr, "        iq = %d row = %d\n", iq, row); }
            rdr[iq] = row;
            if (row != -1)
              { assert((row >= 0) && (row < nrows));
                /* Update the {rChSum} entry of the parent, if exists and is not null: */
                int32_t ip = dt[iq].par;
                if (ip != -1)
                  { assert((ip >= 0) && (ip < iq));
                    rChSum[ip] += (double)row;
                  }
              }
          }
          
        /* Second pass: un-assign and re-assign nodes with no children: */
        for (int32_t iq = ni-1; iq >= 0; iq--) 
          { if (nch[iq] == 0) { unassign_row(iq); } }
        for (int32_t iq = ni-1; iq >= 0; iq--) 
          { if (nch[iq] == 0)
              { assert(rdr[iq] == -1);
                int32_t row = assign_row(iq);
                if (debug) { fprintf(stderr, "        iq = %d row = %d\n", iq, row); }
                rdr[iq] = row;
              }
          }
      }
    /* Find highest used row: */
    int32_t rowMax = -1; /* Max row assigned. */
    for (uint32_t iq = 0;  iq < ni; iq++)
      { if ((rdr[iq] != -1) && (rdr[iq] > rowMax)) { rowMax = rdr[iq]; } }
    /* Reduce plot grid: */
    nrows = rowMax+1;
    fprintf(stderr, "      new cell grid size = %d x %d\n", ncols, nrows);
    free(occ);
    free(nch);
    
    (*ncols_P) = ncols;
    (*nrows_P) = nrows;

    if (debug) { fprintf(stderr, "    < %s\n", __FUNCTION__); }
    return;

    auto int32_t find_free_row(int32_t row_best, int32_t col0, int32_t col1);
      /* Finds the row that has all the cells {col0..col1} free
        and is closest to row {row_best}. */

    int32_t assign_row(int32_t iq)
      { bool_t debug = FALSE;
        if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }
        
        demand(tMin < tMax, "invalid {tMin,tMax}");
        drtree_node_t *q = &(dt[iq]);
        
        if (debug) { fprintf(stderr, "        iq = %d rdr[iq] = %d", iq, rdr[iq]); }
        if (drtree_node_is_null(q))
          { if (debug) { fprintf(stderr, " null\n"); }
            return -1;
          }
        else
          { /* Get its column range: */
            int32_t jloq = dt[iq].tbr - tMin;
            int32_t jhiq = dt[iq].tdt - tMin;
            if (debug) { fprintf(stderr, " cols = {%d .. %d}", jloq, jhiq); }
            assert((0 <= jloq) && (jloq <= jhiq) && (jhiq < ncols));

            /* Choose the ideal row: */
            int32_t rowIdeal;
            int32_t ip = dt[iq].par;
            if ((ip == -1) && (nch[iq] == 0)) 
              { /* No parent or children, place as low as possible: */
                if (debug) { fprintf(stderr, " isolated"); }
                rowIdeal = 0;
              }
            else
              { if (ip != -1)
                  { /* Include the parent row as an attractor: */
                    if (debug) { fprintf(stderr, " ip = %d", ip); }
                  }
                else
                  { /* Include the previously assigned row as atractor: */
                    assert(nch[iq] != 0);
                    if (debug) { fprintf(stderr, " no parent, using ip = iq = %d", iq); }
                    ip = iq;
                  }
                assert((ip >= 0) && (ip <= iq)); /* Sic "<=" */
                drtree_node_t *p = &(dt[ip]);
                assert((0 <= p->tbr) && (p->tbr <= p->tdt) && (p->tdt <= tMax));
                if (debug) { fprintf(stderr, " rdr[ip] = %d", rdr[ip]); }
                assert((rdr[ip] >= 0) && (rdr[ip] < nrows));
                /* Determine the ideal row {rowIdeal} for {iq}, ignoring {occ}: */
                if (debug) { fprintf(stderr, "rChSum[iq] = %.1f nch[iq] = %d", rChSum[iq], nch[iq]); }
                rowIdeal = (int32_t)floor((rChSum[iq] + rdr[ip])/(nch[iq] + 1.0) + 0.5);
              }
            if (debug) { fprintf(stderr, " rowIdeal = %d", rowIdeal); }
            assert((rowIdeal >= 0) && (rowIdeal < nrows));

            /* Choose the row for {iq} that is free and closest to {rowIdeal}: */
            int32_t row = find_free_row(rowIdeal, jloq - 1, jhiq + 1);
            if (debug) { fprintf(stderr, " row = %d\n", row); }
            assert((0 <= row) && (row < nrows));

            /* Mark cells {jloq.(.tdt-tMin)q} on that row as occupied: */
            bool_t *occr = &(occ[row*ncols]);
            for (int32_t col = jloq; col <= jhiq; col++) { occr[col] = TRUE; }

            if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
            return row;
          }
      }

    void unassign_row(int32_t iq)
      {
        int32_t row = rdr[iq];
        if (row == -1) { /* Already unassigned: */ return; }
        assert((row >= 0) && (row < nrows));
        int32_t jloq = dt[iq].tbr-tMin;
        int32_t jhiq = dt[iq].tdt-tMin;
        assert((0 <= jloq) && (jloq <= jhiq) && (jhiq < ncols));
        bool_t *occr = &(occ[row*ncols]);
        for (int32_t col = jloq; col <= jhiq; col++) { occr[col] = FALSE; }
        rdr[iq] = -1;
      }

    int32_t find_free_row(int32_t rowIdeal, int32_t jlo, int32_t jhi)
      { if (jlo < 0) { jlo = 0; }
        if (jhi >= ncols) { jhi = ncols-1; }
        int32_t row = rowIdeal;
        int32_t kd = 0; /* Index into next row increment in {drow}. */
        for (uint32_t kr = 0;  kr < nrows; kr++)
          { /* Try on {row}: */
            bool_t *occr = &(occ[row*ncols]);
            bool_t ok = TRUE;
            for (int32_t col = jlo; (col <= jhi) && ok; col++)
              { if (occr[col]) { ok = FALSE; } }
            if (ok) { return row; }
            
            /* Get to the next closest valid row: */
            row = -1;
            while ((kd < nd) && ((row < 0) || (row >= nrows))) 
              { row = rowIdeal + drow[kd]; kd++; }
            
            /* Must succeed eventually: */
            assert((row >= 0) && (row < nrows));
          }
        demand(FALSE, "failed to find a free row");
      }
  }    
