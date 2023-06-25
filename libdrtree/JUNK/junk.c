
int32_t drtree_add_node(drtree_node_t *dt, int32_t tbr, int32_t tdt, int32_t par);
  /* Adds an individual to the structure {dt}.  Returns
    a sequential index of the individual, starting from 0.
    
    The parameters {tbr} and {tdt} are the birth and death times
    of the individual. Either or both may be negative, but {tdt}
    must not be less than {tbr}.
    
    The parameter {par} is the index of the parent individual, which must have been 
    added previously; or {-1} if the individual is a root. If not {-1},
    then {tbr} must strictly greater than the parent's time of birth,
    and must not exceeds its time of death. */
    
  
drtree_data_t *drtree_data_new(void)
  { 
    drtree_data_t *dt = (drtree_data_t*)notnull(malloc(sizeof(drtree_data_t)), "no mem");
    dt->ni = 0;
    dt->tbr = int32_vec_new(1000);
    dt->tdt = int32_vec_new(1000);
    dt->par = int32_vec_new(1000);
    return dt;
  }
    
void drtree_data_free(drtree_data_t *dt)
  { free(dt->tbr.e);
    free(dt->tdt.e);
    free(dt->par.e);
    free(dt);
  }
  
int32_t drtree_add_node(drtree_data_t *dt, int32_t tbr, int32_t tdt, int32_t par)
  { demand((par == -1) || ((par >= 0) && (par < dt->ni)), "invalid parent index");
    if (par != -1) { demand(tbr > dt->tbr.e[par], "child must be born after parent"); }
    int32_t i = dt->ni;
    int32_vec_expand(&(dt->tbr), i); dt->tbr.e[i] = tbr;
    int32_vec_expand(&(dt->tdt), i); dt->tdt.e[i] = tdt;
    int32_vec_expand(&(dt->par), i); dt->par.e[i] = par;
    dt->ni++;
    return i;
  }

/* INTERNAL STRUCTURES */

typedef struct drtree_compact_info_t 
  { uint32_t nch;      /* Number of children. */
    uint32_t *chi;     /* List of children. */
  } drtree_compact_info_t;
  /* Internal representation and word data
    for an individual in an asexual evolutionary history tree.
    
    The children of the individual are {p.chi[0..p.nch-1]}.
    Any number of children may have the same birth time. */

typedef struct drtree_vis_info_t 
  { /* Parenting info: */
    int32_t jlc;         /* First column of trace. */
    int32_t jhc;         /* last column of trace. */
    int32_t nrc;         /* Number of relevant children. */
    int32_t rpa;         /* Relevant parent index, or {-1}. */
  } drtree_vis_info_t;
  /* Visibility information for an individual {q}, relative to
    a specific plot time range {tMin..tMax}.
  
    The fields {q.jlc} and {q.jhc} are the first and last column of the 
    plot cell grid used by the trace of {q}.  The node {q} is invisible
    if and only if {q.jlc < q.jhc}.

    The field {q.nrc} is the count of relevant children of the individual {q},
    that is, the children whose birth is visible.  In particular, if
    {q} is invisible, then {q.nrc} is zero (but the reciprocal is not true).
    
    The field {q.rpa} is the index of the relevant parent {p} of {q}; namely,
    the parent {p=q.par} of {p} if it exists and if the birth time of {q}
    is in the range {tMin..tMax}.  Otherwise {q.rpa} is {-1}. */

drtree_vis_info_t *drtree_collect_vis_info
  ( int32_t ni,
    drtree_node_t dt[],
    int32_t tMin,
    int32_t tMax
  );
  /* Given an array {dt[0..ni-1]} of basic information for individuals {0..ni-1}, and 
    a time range {tMin..tMax}, returns an array {vf[0..ni-1]} such that 
    {vf[iq]} is the basic visual information about individual {dt[iq]}.
    
    The indices of individuals must have been assigned in topological parenting order,
    meaning that, for each {iq} in {0..ni-1}, {dt[iq].par} is either {-1} or 
    an integer in {0..iq-1}. */

drtree_vis_info_t *drtree_collect_vis_info
  ( int32_t ni,
    drtree_node_t dt[],
    int32_t tMin,
    int32_t tMax
  )
  { 
    drtree_vis_info_t *vf = (drtree_vis_info_t *)notnull(malloc(ni*sizeof(drtree_vis_info_t)), "no mem");
    
    auto void get_node_column_range(drtree_node_t *qd, int32_t *jbr_P, int32_t *jdt_P); 
      /* Given the data record {qd} and visual info record {qv} of an
        individual {q}, computes the range {q.jlc..q.jhc} of cell grid
        columns used by its trace, clipped to {0..ncols-1} where
        {ncols=tMax-tMin+1}. Returns the bounds in {*jbr_P} and {*jdt_P}.
        Note that the range may be empty. */

    /* Scan in chrono order, initializing {p.nco} to the life span of the node only: */
    for (int32_t iq = 0; iq < ni; iq++) 
      { drtree_vis_info_t *vfq = &(vf[iq]);
        vfq->nrc = 0;
        vfq->rpa = -1;
        get_node_column_range(&(dt[iq]), &(vfq->jlc), &(vfq->jhc));
      }
      
    /* Scan in reverse chrono order, defining the {nrc} and {rpa} fields: */
    for (int32_t iq = ni-1; iq >= 0; iq--)
      { drtree_node_t *dtq = &(dt[iq]);
        drtree_vis_info_t *vfq = &(vf[iq]);
        int32_t ip = dtq->par;
        if ((ip != -1) && (dtq->tbr >= tMin) && (dtq->tbr <= tMax)) 
          { /* Node {iq} is a relevant child of {ip}. */
            assert((ip >= 0) && (ip < iq));
            drtree_node_t *dtp = &(dt[ip]);
            drtree_vis_info_t *vfp = &(vf[ip]);
            assert(dtp->tbr < dtq->tbr); /* Can't be parent at or before birth. */
            assert(dtp->tdt >= dtq->tbr); /* Can't be parent after death. */
            vfp->nrc++;
            vfq->rpa = ip;
          }
      }

    return vf;

    void get_node_column_range(drtree_node_t *dtq, int32_t *jlc_P, int32_t *jhc_P)
      {
        demand(tMin <= tMax, "invalid time range");
        int32_t ncols = tMax - tMin + 1;
        assert(dtq->tbr <= dtq->tdt);
        int32_t jlc = dtq->tbr - tMin; if (jlc < 0) { jlc = 0; }
        int32_t jhc = dtq->tdt - tMin; if (jhc >= ncols) { jhc = ncols - 1; }
        (*jlc_P) = jlc;
        (*jhc_P) = jhc;
      }
  }

void drtree_check_vis_info
  ( int32_t ni,
    drtree_node_t dt[],
    int32_t tMin,
    int32_t tMax,
    drtree_vis_info_t vf[]
  );
  /* Checks the consistency of {vf[0..ni-1]} given {dt[0..ni-1]},
    {tMin}, and {tMax}. */

void drtree_check_vis_info
  ( int32_t ni,
    drtree_node_t dt[],
    int32_t tMin,
    int32_t tMax,
    drtree_vis_info_t vf[]
  )
  {
    int32_t ncols = tMax - tMin + 1;
    
    /* Check {jlc,jhc}, and {rpa,nrc} of invisible nodes: */
    for (int32_t iq = 0; iq < ni; iq++)
      { drtree_node_t *dtq = &(dt[iq]);
        assert(dtq->tbr <= dtq->tdt);
        drtree_vis_info_t *vfq = &(vf[iq]);
        assert(vfq->jlc == (int32_t)imax(0, dtq->tbr - tMin));
        assert(vfq->jhc == (int32_t)imin(ncols-1, dtq->tdt - tMin));
        if (vfq->jlc > vfq->jhc)
          { /* Node {iq} is invisible: */
            assert(vfq->rpa == -1);
            assert(vfq->nrc == 0);
          }
      }
      
    /* Check {nrc,rpa}: */
    int32_t nrc[ni];
    for (int32_t iq = 0; iq < ni; iq++) { nrc[iq] = 0; }
    for (int32_t iq = ni-1; iq >= 0; iq--)
      { drtree_node_t *dtq = &(dt[iq]);
        drtree_vis_info_t *vfq = &(vf[iq]);
        if (vfq->jlc <= vfq->jhc)
          { /* {q} is partly visible */
            assert(vfq->nrc == nrc[iq]);
            if ((dt[iq].tbr >= tMin) && (dt[iq].tbr <= tMax))
              { /* Birth of {q} is visible. */
                int32_t ip = dtq->par;
                if (ip == -1)
                  { assert(vfq->rpa == -1); }
                else
                  { assert((ip >= 0) && (ip < iq));
                    drtree_node_t *dtp = &(dt[ip]);
                    drtree_vis_info_t *vfp = &(vf[ip]);
                    assert((dtp->tbr < dtq->tbr) && (dtp->tdt >= dtq->tbr));
                    assert(vfp->jlc <= vfp->jhc);
                    assert(vfq->rpa == ip);
                    nrc[ip]++;
                  }
              }
            else
              { /* Birth of {q} is invisible: */
                assert(dtq->tbr < tMin);
                assert(vfq->rpa == -1);
              }
          }
        else
          { /* {q} is invisible: */
            assert(vfq->nrc == 0);
            assert(vfq->rpa == -1);
          }
      }
  }

void drtree_get_subtree_column_range
  ( drtree_node_t *qd,
    int32_t tMin,
    int32_t tMax,
    int32_t *jlo_P,
    int32_t *jhi_P
  );
  /* Given the data record {qd} and visual info record {qv} of an
    individual {q}, computes the range {jlo..jhi} of cell grid columns
    used by the traces of {q} and all its relevant descendants, clipped
    to {0..ncols-1} where {ncols=tMax-tMin+1}. Returns the bounds in
    {*jlo_P} and {*jhi_P}. 
    
    Note that the range may be empty. In particular, if {q} is
    invisible, the range will be empty, even if some descendants
    are visible -- since they will not be /relevant/ descendants. */
 
void drtree_get_subtree_column_range
  ( drtree_node_t *qd,
    int32_t tMin,
    int32_t tMax,
    int32_t *jlo_P,
    int32_t *jhi_P
  )
  {
    demand(tMin <= tMax, "invalid time range");
    if (qv->nco == 0)
      { /* Node {q} is invisible hence irrelevant: */
        assert((qd->tdt < tMin) || (qd->tbr >= tMax));
        (*jlo_P) = 1; (*jhi_P) = 0;
      }
    else
      { /* Node {q} is at least partly visible: */
        /* Get first visible column {jbr} of {q}: */
        int32_t ncols = tMax - tMin + 1;
        int32_t jbr = qd->tbr - tMin; if (jbr < 0) { jbr = 0; }
        /* Visible subtree extends {qv->nco} cols from {jbr}: */
        int32_t jlo = jbr;
        int32_t jhi = jlo + qv->nco - 1;
        assert(jhi < ncols);
        (*jlo_P) = jlo;
        (*jhi_P) = jhi;
      }
  }
   
drtree_compact_info_t *drtree_compact_collect_children(drtree_data_t *dt)
  { 
    int32_t ni = dt->ni;
    /* Collectr children counts: */
    drtree_compact_info_t *ci = (drtree_compact_info_t *)notnull(malloc(ni*sizeof(drtree_compact_info_t)), "no mem");
    for (int32_t iq = 0; iq < ni; iq++) 
      { ci[iq].nch = 0; ci[iq].chi = NULL;
        int32_t p = dt->{par|rpa}.e[iq];
        if (p != -1) { ci[p].nch++; }
      }
    /* Alocate the children lists: */
    for (int32_t iq = 0; iq < ni; iq++) 
      { ci[iq].chi = (int32_t*)notnull(malloc(ci[iq].nch*sizeof(int32_t)), "no mem");
        ci[iq].nch = 0;
      }
    /* Collect the children: */
    for (int32_t iq = 0; iq < ni; iq++) 
      { int32_t p = dt->{par|rpa}.e[iq];
        if (p != -1) 
          { ci[p].chi[ci[p].nch] = iq;
            ci[p].nch++; 
          }
      }
    return ci;
  }


void foo(void)
  { 
    /* Initialize {rdr[iq]} to {-1} for all {iq}.  That will mean "not assigned yet". */
    for (int32_t iq = 0; iq < ni; iq++) { rdr[iq] = -1; }
    
    auto int32_t assign_lineage_rows(int32_t iq);
      /* Sets {rdr[iq]}, and {rdr[jq]} for all descendants of {iq}.
        Updates the horizon {hir[0..ncols-1]} as needed. 
        Returns the number of individuals assigned.*/
    
    /* Scan individuals in chrono order, assign lineages of unassigned ones: */
    int32_t rowMax = -1;
    for (int32_t iq = 0; iq < ni; iq++)
      { if (rdr[iq] == -1)
          { /* Still unassigned: */
            assert(st->{par|rpa}.e[iq] == -1); /* Otherwise should have been assigned. */
            int32_t nass = assign_lineage_rows(iq);
            if (debug) { fprintf(stderr, "        root indiv %d assigned to row %d, %d descendants\n", iq, rdr[iq], nass); }
            assert(rdr[iq] >= 0);
          }
        if (rdr[iq] > rowMax) { rowMax = rdr[iq]; }
      }
    
    /* Reduce plot grid: */
    nrows = rowMax+1;
    fprintf(stderr, "      new cell grid size = %d x %d\n", ncols, nrows);

    free(nch);
    (*nrows_P) = nrows;
    if (debug) { fprintf(stderr, "    < %s\n", __FUNCTION__); }
    return;
    
    auto int32_t assign_to_lowest_free_row(int32_t col0, int32_t col1);
      /* Finds the lowest row where it is possible to place a 
        a trace that spans grid columns {col0...col1} so that it is 
        above all previous traces that enter those columns
        and has at least one free cell on each side.
        
        Also updates {hir[col0..col1]} to that row. */
        
    auto bool_t child_goes_below(int32_t r, int32_t nchi);
      /* Decides whether child {r} out of {nchi} children 
        should go below or above the parent: */
    
    int32_t assign_lineage_rows(int32_t iq)
      { bool_t debug = FALSE;
        if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }
    
        /* Get main parameters of individual {iq}: */
        int32_t tbri = st->tbr.e[iq];
        int32_t tdti = st->tdt.e[iq];
        int32_t nchi = st->nch.e[iq];
        if (nchi == -1) { assert(tdti == -1); nchi = 0; tdti = tbri; }
        
        /* Determine the column range of {iq}'s trace, including the time of death: */
        int32_t col0 = tbri - tMin;
        int32_t col1 = tdti - tMin;
        
        /* Assign its lineage: */
        int32_t nass = 1; /* Number of individuals in lineage. */
        int32_t rowi = -1;
        if (nchi == 0)
          { /* Just place {iq} as low as possible. */
            rowi = assign_to_lowest_free_row(col0, col1);
          }
        else
          { /* Find all the chidlren. We need not scan too far since life is short. */
            int32_t chi[nchi];
            int32_t kchi = 0; /* Counts children found. */
            for (int32_t jq = 0; (jq < ni) && (kchi < nchi); jq++)
              { int32_t p = st->{par|rpa}.e[jq];
                if (p == iq) { chi[kchi] = jq; kchi++; }
              }
            assert(kchi == nchi);
             
            /* Assign some children in order of birth: */
            for (int32_t r = 0; r < nchi; r++)
              { if (child_goes_below(r, nchi))
                  { int32_t nassj = assign_lineage_rows(chi[r]);
                    nass += nassj;
                  }
              }

            /* Tentatively assign {iq} as low as possible: */
            rowi = assign_to_lowest_free_row(col0, col1);
            
            /* Assign the remaining children in reverse order: */
            for (int32_t r = nchi-1; r >= 0; r--)
              { if (! child_goes_below(r, nchi))
                  { int32_t nassj = assign_lineage_rows(chi[r]);
                    nass += nassj;
                  }
              }
              
            /* Raise {iq} to midway of its children: */
            fprintf(stderr, "!! individual %d -- raising not implemented\n", iq);
          }
        rdr[iq] = rowi;
        if (debug) { fprintf(stderr, "          indiv %d assigned to row %d, %d children, %d descendants \n", iq, rowi, nchi, nass); }
        if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
        return nass;
      }
      
    bool_t child_goes_below(int32_t r, int32_t nchi)
      { /* Let's try to put 1, even below: */
        return (r == 1) || ((r % 2) == 0);
      }
    
    int32_t assign_to_lowest_free_row(int32_t col0, int32_t col1)
      { assert((0 <= col0) && (col0 <= col1) && (col1 < ncols));
        /* Find maximum of {hir[col0-1..col1+1} to leave space on each side: */
        int32_t hirMax = -1;
        for (int32_t col = col0-1; col <= col1+1; col++)
          { if ((col >= 0) && (col < ncols))
              { if (hir[col] > hirMax) { hirMax = hir[col]; } }
          }

        /* Assign to row above: */
        int32_t row = hirMax + 1;
        assert((row >= 0) && (row < nrows));

        /* Update horizon on cols {col0..col1}: */
        for (int32_t col = col0; col <= col1; col++) { hir[col] = row; }
        
        return row;
      }
  }

    
    /* The top horizon table: */
    int32_t hir[ncols];
    for (int32_t col = 0; col < ncols; col++) { hir[col] = -1; }
