/* See {drtree_planar.h} */
/* Last edited on 2024-11-04 07:31:24 by stolfi */

#define drtree_planar_C_COPYRIGHT \
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
#include <jsrandom.h>
#include <in.h>
#include <vec.h> 

#include <drtree.h>
#include <drtree_planar.h>

/* INTERNAL STRUCTURES */

typedef struct drtree_planar_info_t 
  { int32_t nch;         /* Number of children. */
    int32_t *chi;        /* List of indices of children. */
    /* Subtree shape info: */
    int32_t jlo;         /* First visible column of subtree. */
    int32_t jhi;         /* Last visible column of subtree. */
    int32_t *uph;        /* Upper horizon of subtree including parent. */
    int32_t upMax;       /* Max of {uph[0..jhi-jlo} */
    int32_t *loh;        /* Lower horizon of subtree. */
    int32_t loMax;       /* Max of {loh[0..jhi-jlo]}. */
    /* Placement of subtree relative to parent: */
    int32_t row;         /* Relative assigned row of the plot, or 0. */
    bool_t flp;          /* True if the subtree is to be flipped. */
  } drtree_planar_info_t;
  /* Internal complementary work data for some individual {iq} in an 
    asexual evolutionary history tree.
    
    The following comments assume that {q} is a non-null node, that
    {q.tbr-tMin .. q.tdt-tMin} is contained in {0..ncols}, and that
    {q.nch} is the number of children of {q}.

    CHILDREN LIST
    
    If {nch[iq]} is zero, {q.chi} is {NULL}; otherwise {q.chi[0..nch[iq]-1]} are the
    indices of those relevant children.
    
    COLUMNS SPANNED BY SUBTREE

    The fields {q.jlo} and {q.jhi} are the first and last cell grid
    columns spanned by the subree rooted at {q}. That is, {q.jlo} is
    {q.tbr-tMin}, and {q.jhi} is is the max value of {d.tdt-tMin} for
    any descendant {d} of {q}, including {q} itself.
    
    In the comments that follows, {q.} means either {af[iq].} or {dt[iq].}, where 
    {iq} is the ID of node {q}.  Also {q.nco} is a shorthand of {q.jhi-q.jlo+1}, the
    number of columns spanned by the subtree rooted at {q}. 
     
    If {q} has a parent {p = q.par}, then {q.row} is the row which {q}
    was assigned to, but relative to {p}'s assigned row. This number may
    be positive or negative but not zero. If {q} is a root, its parent's
    row is assumed to be {-1} for this purpose; so {q.row} is the
    absolute row assigned to {q}, plus one -- also never zero. While not
    assigned, {q.row} is set to zero.
    
    The field {q.uph}, if not null, is the /upper horizon/ of the
    drawing of the tree rooted at {q}. It is an array of {q.nco}
    integers such that {q.uph[k]} is the highest row ocupied by the
    subtree of {q} on column {q.jlo+k}. These horizon values are
    relative to {q}'s assigned row, with positive values meaning
    `above'. The field {q.upMax} is the maximum of {q.uph[0..q.nco-1]}.
    
    The field {q.loh}, if not null, is the /lower horizon/, defined similarly
    to {q.uph} except that it tells the lowest occupied row, and positive 
    values mean `below'.  The field {q.loMax} is likewise the maximum 
    of {q.loh[0..q.nco-1]}.
    
    The field {q.flp}, if true, says that the whole subtree rooted at
    the individual {q} should be drawn flipped upside-down around the
    individual itself. That is, the the upper and lower horizons
    {q.uph,q.loh} sould be swapped, like {q.upMax} and {q.loMax}, each
    child {c = q.chi[k]} should be drawn on row {p.row-c.row}, rather
    than {p.row+c.row}, and its bit {c.flp} shoud be considered inverted
    as well.
    
    Thus, if {q.flp} is false (taking into account the {flp} fields of
    {q}'s ancestors), the subtree rooted at {q} is contained in the
    rectangle spanning columns {q.jlo .. q.jhi} and rows {i-q.loMax
    .. i+q.upMax} where {i} is the actual row assigned to {q}. More
    precisely, on each column {q.jlo+k}, with {k} in {0..q.nco-1}, the
    subtree is contained in the range {i-q.loh[k] .. i+q.uph[k]}.
    If {q.flp} is true (taking into account the {flp} fields of
    {q}'s ancestors) the same is true with {uph,loh} swapped and
    {upMax,loMax} swapped. */
  
#define drtree_planar_NO_ROW (- drtree_indivs_MAX) 
  /* A negative number less than any valid relative or absolute row index.
    Cannot de {INT32_MIN} because negation would overflow. */
     
#define drtree_planar_enum_sub_num_MAX 8
  /* Max number {nch} of subtrees for exhaustive enumeration. Must be 
    such that {2^{2*nch-1}} will not overflow. */
    
#define drtree_planar_enum_sub_try_MAX 32
  /* Max number of above/below ad flip arrangements to enumerate 
    exhaustively when placing the subtrees of a tree. If {2^{2*nch-1}}
    is more than this, will do a partial enumeration. */
       
#define drtree_planar_enum_top_num_MAX 3
    /* Max number {nr} of top level trees for exhaustive enumeration. Must be 
    such that {nr!*2^{nr-1}} will not overflow. */
   
#define drtree_planar_enum_top_try_MAX 24
  /* Max number of of permutation ad flip arrangements to enumerate exhaustively when placing the 
    top-level trees.  If {nr!*2^{nr-1}} is greater than this, does . */

/* INTERNAL PROTOTYPES */

drtree_planar_info_t *drtree_planar_info_collect
  ( int32_t ni, 
    drtree_node_t dt[],
    int32_t nch[],
    int32_t tMin, 
    int32_t tMax
  );
  /* Returns an array {af[0..ni-1]} of {drtree_planar_info_t} 
    based on the data {dt[0..ni-1]}.  The fields {af[iq].jlo)}, {af[iq].jhi}
    and {af[iq].nch} will be set, all other fields will be undefined.
    
    The life span {q.tbr..q.tdt} of every node {q=dt[iq]} must be a
    subset of {tMin..tMax}. Assumes that {nch[iq]} is the child count of
    node {q}. */

void drtree_planar_info_free(int32_t ni, drtree_planar_info_t *af);
  /* Releases the storage used by the record {af[0..ni-1]} and all its internal storage,
    including horizons and child lists. */

void drtree_planar_sort_children(int32_t ni, drtree_node_t dt[], int32_t nch, int32_t chi[]);
  /* Sorts the list {chi[0..nch-1]} of individual IDs in decreasing order of their
    birth times, as specified in the table {dt[0..ni-1]}. */

void drtree_planar_place_subtrees_relative
  ( int32_t iq,
    int32_t ni,
    drtree_node_t dt[],
    drtree_planar_info_t af[]
  );
  /* Given data tables {dt[0..ni-1]} and {af[0..ni-1]}, chooses a
    relative placement for the {q.nch} subtrees of the node {q = dt[iq]},
    by setting the fields {c.flp} and {c.row} of each child
    {c = dt[ic]}, where {ic = q.chi[r]} for {r} in {0..nch-1}.
    
    Assumes that the children are sorted in increasing order of birth
    time, and that the horizons {c.uph[0..c.nco-1]} and
    {c.loh[0..c.nco-1]} of each child are defined.
    
    Also defines the upper and lower horizons {af[iq].uph[0..q.nco-1]}
    and {af[iq].loh[0..q.nco-1]}. */

int32_t drtree_planar_place_top_trees
  ( int32_t ni,
    drtree_node_t dt[],
    drtree_planar_info_t af[],
    int32_t nr,
    int32_t rut[],
    int32_t tMin,
    int32_t ncols,
    int32_t rdr[]
  );
  /* Given the list {rut[0..nr-1]} of the 
    indices of the roots of the evolution forest, chooses an
    absolute row index {rdr[iq]} for each root {iq = rut[r]}, so as to
    avoid overlap of the trees.  Returns the number of rows
    used by the diagram. 
    
    Assumes defined the horizons {af[iq].uph[0..q.nco-1]} and
    {af[iq].loh[0..q.nco-1]} of the tree rooted at each root {iq}.  */

void drtree_planar_set_absolute_rows
  ( int32_t ni,
    drtree_node_t dt[],
    drtree_planar_info_t af[],
    int32_t rdr[]
  );
  /* For {iq} in {0..ni-1}m sets {rdr[iq]} to the row of the  diagram assigned to the trace of 
    node {q = dt[iq]}.  Also propagates the {flp} bits so that they are absolute rather than 
    relative to the parent.
    
    If {q = dt[iq]} is a root node,  assumes that {af[iq].row} is the absolute 
    diagram row. If {q} has a parent {p}, assumes that {q.row} is relative to the 
    absolute row of {p}, either up or down depending on the subtree of {p} is flipped
    relative to the world. If {q} is a null node, sets {rdr[iq]} to {-1}. */

typedef void drtree_planar_abo_flp_proc_t(bool_t abo[], bool_t flp[]);
  /* Type of a procedure that processes one combination of
    {abo[0..nch-1]} and {flp[0..nch-1]} generated by 
    {drtree_planar_enum_subtree_placements}. */

void drtree_planar_enum_subtree_placements
  ( int32_t nch,
    drtree_planar_abo_flp_proc_t try_abo_flp
  );
  /* Enumerates possible placements of the subtrees of {nch} children of an individual
    relative to the parent.  
    
    If {nch} is small enough, enumerates all {2^{2*nch-1}} combinations
    of above;below and flipping flags {abo[0..nch-1]} and
    {flp[0..nch-1]} and calls {try_abo_flp(abo,flp)} for each of them.
    The element {flp[nch-1]} will always be false, since negating all
    {abo} and {flp} flags is essentially the same arrangement.
    
    If {nch} is too big, enumerates only a few combinations, at random,
    starting with all even children below, all odd-children above. */

typedef void drtree_planar_ord_flp_proc_t(int32_t ord[], bool_t flp[]);
  /* Type of a procedure that processes one combination of
    {ord[0..nr-1]} and {flp[0..nr-1]} generated by 
    {drtree_planar_enum_top_tree_placements}. */

void drtree_planar_enum_top_tree_placements
  ( int32_t nr,
    drtree_planar_ord_flp_proc_t try_ord_flp
  );
  /* If {nr} is small enough, enumerates all {nr!*2^{nr-1}} combinations
    of root order {ord[0..nr-1]} and flipping flags 
    {flp[0..nr-1]} and calls {try_ord_flp(ord,flp)} for each of them.
    
    The element {flp[nr-1]} is always false, since reversing 
    the order {ord} of all roots and negating all {flp}
    bits should yield an arrangement with the same height.
    
    If {nr} is too big, enumerates only a few combinations, at random,
    starting with bigger trees first. */

void drtree_planar_check_info
  ( int32_t ni,
    drtree_node_t dt[],
    drtree_planar_info_t af[],
    int32_t tMin,
    int32_t tMax
  );
  /* Checks the consistency of {af[0..ni-1]} given {dt[0..ni-1]},
    {tMin}, and {tMax}.  Should be called after {nch,ch} have been
    defined but before nodes have been placed (i.e. before {uph,loh,row,flp}
    have been computed).  Assumes that {dt[0..ni-1]} is OK (see {drtree_check_nodes}).  */

/* IMPLEMENTATIONS */ 

void drtree_planar_arrange
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

    drtree_check_nodes(ni, dt, tMin, tMax);
    
    /* Get the child counts: */
    int32_t *nch = drtree_count_children(ni, dt);
    
    int32_t ncols = tMax - tMin + 1;
    int32_t nrows = ni; /* Tentative -- to be reduced later. */
    
    /* The general strategy is to recursively assign rows to the lineage
      of each individual.  
      
      Using the fields {q.tbr-tMin}, {q.tdt-tMin}, {nch[iq]}, and {q.par} from
      {af[0..ni-1]}, the procedure first collects the relevant children
      {q.chi} are first set for all nodes, by a sweep in chronological
      order (parent before children).
      
      Nodes are then processed in reverse chronological order (children
      before parent). Once each node {q} is processed, its subtree shape
      fields {q.uph,q.upMax,q.loh,q.loMax} are defined and fixed.
      When its parent {p} is processed, the parent-relative layour
      fields {q.row} and {q.flp} are defined and fixed.
      
      When processing a visible node {q}, the procedure decides for
      each relevant child {c} with a visible subtree (that is, with {(c.tbr-tMin) <= (c.tdt-tMin)})
      whether it should be placed above or below {c} 
      and whether it should be flipped. It also decided in which order the
      subtrees on each side of {q} should be stacked.
      
      In the worst case, there are {(m+1)2^m m!} possibilities for these
      choices, where {m} is the number of visible subtrees. Some of
      those choices are ruled out by the planarity requirement, and some
      can be eliminated because they have no effect on the resulting
      layout. Anyway, if the number of possibilities is to large, the
      procedure will use some heuristic and try only a subset of all
      possible choices.
      
      For each of those possible combinations of flipping, above/below,
      and ordering, the visible child subtrees that are chosen to go
      above {q} will be stacked (possibly flipped) above the trace of
      {q}, and those chosen to go below will be stacked (possibly
      flipped) below that trace, as snugly as possible, without
      overlaps.
       
      Actually, if {m>0}, half of the possibilities can be ignored,
      because the whole arrangement can later be flipped
      around the parent {q}. Also, if {c} is visible but nas no visible 
      children, flipping its subtree has no effect, so it contributes
      only two possibilities (above or below {q}) rather than 4.
     
      Moreover, the subtrees of the children {c} must be stacked first,
      in order of decreasing birth time, to ensure planarity.
      
      For each combination of flipping, above/below, and ordering that
      is considered, the procedure will compute the resulting horizons
      {q.uph} and {q.loh}, and the total number of cells occupied by the
      subtree of {q}, which is {SUM{ uph[k]+loh[k]+1 : k \in
      {0..q.jhi-q.jlo}}}. Of all those possibilities that are considered,
      the procedure will choose the one which results in the minimum
      radius. */
      
    /* Allocate the info records: */
    drtree_planar_info_t *af = drtree_planar_info_collect(ni, dt, nch, tMin, tMax);
    drtree_planar_check_info(ni, dt, af, tMin, tMax);
          
    /* Now scan in reverse chrono order and define the placement 
      of subtrees relative to parent. Also count the roots: */
    int nr = 0; /* Count of roots. */
    for (int32_t iq = ni-1; iq >= 0; iq--)
      { drtree_node_t *q = &(dt[iq]);
        drtree_planar_info_t *afq = &(af[iq]);
        if (drtree_node_is_null(q))
          { /* Ignore node: */
            assert(afq->chi == NULL);
            assert(afq->jhi < afq->jlo);
            assert(nch[iq] == 0);
          }
        else      
          { /* Choose placement of subtrees relative to {q}: */
            drtree_planar_place_subtrees_relative(iq, ni, dt, af);
            assert((afq->uph != NULL) && (afq->loh != NULL)); 
            if (q->par == -1) { nr++; }
          }
      }
      
    /* Now collect all roots rut[0..nr-1]: */
    int32_t rut[nr];
    { int32_t kr = 0; /* Roos found so far. */
      for (int32_t iq = 0; iq < ni; iq++)
        { drtree_node_t *q = &(dt[iq]);
          if ((! drtree_node_is_null(q)) && (q->par == -1))
            { rut[kr] = iq; kr++;
              drtree_planar_info_t *afq = &(af[iq]);
              assert((afq->uph != NULL) && (afq->loh != NULL)); 
            }
        }
      assert(kr == nr);
    }
    
    /* Now place all roots absolutely, get the diagram's height: */
    nrows = drtree_planar_place_top_trees(ni, dt, af, nr, rut, tMin, ncols, rdr);
    
    /* Now propagate the placement to the descendants: */
    drtree_planar_set_absolute_rows(ni, dt, af, rdr);
    
    /* Cleanup: */
    drtree_planar_info_free(ni, af);
    
    (*ncols_P) = ncols;
    (*nrows_P) = nrows;
    
    return;    
  }
 
drtree_planar_info_t *drtree_planar_info_collect
  ( int32_t ni, 
    drtree_node_t dt[],
    int32_t nch[], 
    int32_t tMin, 
    int32_t tMax
  )
  { bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }

    assert((ni >= 0) && (ni < drtree_indivs_MAX));
    assert(tMin < tMax);

    drtree_planar_info_t *af = (drtree_planar_info_t*)notnull(malloc(ni*sizeof(drtree_planar_info_t)), "no mem"); 
    
    /* Scan in chrono order, setting {q.nch}, allocating {q.chi}, and 
      initializing {q.jhi} and {q.jlo} to the life span of the node only: */
    for (int32_t iq = 0; iq < ni; iq++) 
      { if (debug) { fprintf(stderr, "        iq = %d", iq); }
        drtree_node_t *q = &(dt[iq]);
        drtree_planar_info_t *afq = &(af[iq]);
        if (! drtree_node_is_null(q))
          { afq->nch = nch[iq];
            afq->chi = (afq->nch == 0 ? NULL : (int32_t *)notnull(malloc(nch[iq]*sizeof(int32_t)), "no mem"));
            afq->jlo = q->tbr - tMin;
            afq->jhi = q->tdt - tMin;
            if (debug) { fprintf(stderr, " nch = %d span {%d..%d} -> {%d..%d}", afq->nch, q->tbr, q->tdt, afq->jlo, afq->jhi); }
           }
        else
          { assert(nch[iq] == 0);
            assert(q->par == -1);
            afq->nch = 0;
            afq->chi = NULL;
            afq->jlo = 1;
            afq->jhi = 0;
            if (debug) { fprintf(stderr, " null"); }
          }
          
        afq->uph = NULL;
        afq->upMax = drtree_planar_NO_ROW;
        afq->loh = NULL;
        afq->loMax = drtree_planar_NO_ROW;
        afq->row = drtree_planar_NO_ROW;
        afq->flp = FALSE;
        if (debug) { fprintf(stderr, "\n"); }
        
        nch[iq] = 0; /* Temporarily, to be restored. */
      }

    /* Scan in reverse chrono order, including the column extent
      of children in that of the  parent and filling its {chi} list: */
    if (debug) { fprintf(stderr, "\n"); }
    for (int32_t iq = ni-1; iq >= 0; iq--)
      { if (debug) { fprintf(stderr, "        iq = %d", iq); }
        drtree_node_t *q = &(dt[iq]);
        if (! drtree_node_is_null(q))
          { if (debug) { fprintf(stderr, " tbr = %d", q->tbr); }
            drtree_planar_info_t *afq = &(af[iq]);
            assert((afq->jhi >= afq->jlo) && (q->tbr >= tMin) && (q->tdt <= tMax));
            int32_t ip = q->par;
            if (debug) { fprintf(stderr, " parent = %d", ip); }
            if (ip != -1)
              { assert((ip >= 0) && (ip < iq));
                drtree_node_t *p = &(dt[ip]);
                if (debug) { fprintf(stderr, " span {%d..%d}", p->tbr, p->tdt); }
                assert(! drtree_node_is_null(p)); /* A null node can't be parent. */
                assert(p->tbr < q->tbr); /* Can't be parent at or before birth. */
                assert(p->tdt >= q->tbr); /* Can't be parent after death. */
                drtree_planar_info_t *afp = &(af[ip]);
                assert(afp->jhi >= afp->jlo);
                assert((afq->jlo >= afp->jlo) && (afq->jlo <= afp->jhi));

                /* Update {p.jlo..p.jhi} to include the subtree of {iq}: */
                afp->jhi = (int32_t)imax(afp->jhi, afq->jhi);
                if (debug) { fprintf(stderr, " col span expanded to {%d..%d}", afp->jlo, afp->jhi); }
                
                /* Add {q} to the children or {p}: */
                assert(afp->chi != NULL);
                afp->chi[nch[ip]] = iq;
                nch[ip]++; /* This should eventually restore the {nch} fields. */
              }
            else
              { if (debug) { fprintf(stderr, " root"); } }
          }
        else
          { if (debug) { fprintf(stderr, " null"); } }
        if (debug) { fprintf(stderr, "\n"); }
      }
      
    /* Now scan in chrono order sorting the children, just in case: */
    for (int32_t iq = 0; iq < ni; iq++)
      { drtree_node_t *q = &(dt[iq]);
        if (! drtree_node_is_null(q))
          { assert ((q->tbr >= tMin) && (q->tbr <= tMax));
            drtree_planar_info_t *afq = &(af[iq]);
            assert (afq->jhi >= afq->jlo);
            assert(nch[iq] == afq->nch);
            drtree_planar_sort_children(ni, dt, afq->nch, afq->chi);
          }
      }
      
    if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
    return af;
  }
 
void drtree_planar_sort_children(int32_t ni, drtree_node_t dt[], int32_t nch, int32_t chi[])
  { 
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }
    
    assert((ni >= 0) && (ni < drtree_indivs_MAX));

    for (int32_t r = 1; r < nch; r++)
      { int32_t ir = chi[r];
        int32_t tbrr = dt[ir].tbr;
        int32_t s = r; 
        while ((s > 0) && (dt[chi[s-1]].tbr < tbrr))
          { chi[s] = chi[s-1]; s--; }
        chi[s] = ir;
      }
    if (debug) 
      { fprintf(stderr, "        chi.tbr = ("); 
        for (int32_t r = 0; r < nch; r++) { fprintf(stderr, " %d", dt[chi[r]].tbr); }
        fprintf(stderr, " )\n"); 
      }
    if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
  }

void drtree_planar_place_subtrees_relative
  ( int32_t iq,
    int32_t ni,
    drtree_node_t dt[],
    drtree_planar_info_t af[]
  )
  { /* bool_t debug = FALSE; */
    bool_t debug = (iq == 1);
    if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }

    assert((ni >= 0) && (ni < drtree_indivs_MAX));

    drtree_node_t *q = &(dt[iq]);
    drtree_planar_info_t *afq = &(af[iq]);
    int32_t nch = afq->nch;
    int32_t ncoq = afq->jhi - afq->jlo + 1;
         
    if (debug) { fprintf(stderr, "        iq = %d nch = %d ncoq = %d\n", iq, nch, ncoq); }

    /* Allocate horizon lists: */
    assert(afq->uph == NULL);
    afq->uph = (int32_t *)notnull(malloc(ncoq*sizeof(int32_t)), "no mem");
    assert(afq->loh == NULL);
    afq->loh = (int32_t *)notnull(malloc(ncoq*sizeof(int32_t)), "no mem");
    
    /* Best above/below and flip arrangement found: */
    bool_t abo_best[nch]; /* Tells which children go above the parent. */          
    bool_t flp_best[nch]; /* Tells which subtrees should be flipped up/down */         
    int32_t radius_best = INT32_MAX; /* Radius extent of {abo_best,flp_best}. */ 
    
    auto void try_abo_flp(bool_t abo[], bool_t flp[]);
      /* To be given to {drtree_planar_enum_subtree_placements}. Assumes
        {abo[0..nch-1]} are above/below for the children of {q}, and
        {flp[0..nch-1]} are their up/down flip assignments. Computes the
        relative row assignments for that arrangement, the resulting
        horizons of {q}, and the radius (maximum deviation
        of the horizons from the parent row); and
        updates {abo_best,flp_best,radius_best} accordingly. */
    
    auto void assign_relative_child_rows(bool_t abo[], bool_t flp[]);
      /* For each child {c = dt[ic]} of {q}, where {ic = q.chi[r]} with {r} in {0..q.nch-1},
        sets the flip bit {c.flp} of of {c} to {flp[r]}. Then computes
        the relative row displacement {af[ic].row} of each child
        relative to {q}, by stacking the subtrees of the children (each
        flipped as indicated) above or below {q}, as specified by
        {abo[r]}. Also defines the horizons {q.uph[0..q.nco-1]} and
        {q.loh[0..q.nco-1]}.
        
        Assumes that the children are sorted in increasing order of
        birth time, and that the horizons {c.uph[0..c.nco-1]} and
        {c.loh[0..c.nco-1]} of each child {c} are defined. */
    
    auto int32_t compute_envelope_radius(void);
      /* Computes the ares of the subtree of {q}, that is,
        the number of grid cells in the interval
        {-q.loh[j] .. +q.uph[j]} for {j} in 
        {0..q.jhi-q.jlo}.  Assumes that these horizons are 
        defined. */

    auto void free_child_horizons(void);
      /* Frees the horizon lists of all children of {q}. */

    /* Enumerate {abo,flp}, find best: */
    drtree_planar_enum_subtree_placements(nch, try_abo_flp);
     
    /* Set best placement, compute final horizons of {q}: */
    assign_relative_child_rows(abo_best, flp_best);
    
    /* Cleanup horizons of children: */
    free_child_horizons();
    
    if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
        
    return;
    
    /* INTERNAL IMPLEMENTATIONS */
    
    void try_abo_flp(bool_t abo[], bool_t flp[])
      { assign_relative_child_rows(abo, flp);
        int32_t radius = compute_envelope_radius();
        if (radius < radius_best)
          { for (int32_t r = 0; r < nch; r++)
              { abo_best[r] = abo[r];
                flp_best[r] = flp[r];
                radius_best = radius; 
              }
          }
      }   
      
    void assign_relative_child_rows(bool_t abo[], bool_t flp[])
      { /* Initialize {q}'s horizons with {q}'s trace: */
        if (debug) { fprintf(stderr, "          > %s\n", __FUNCTION__); }
        if (debug) { fprintf(stderr, "            iq = %d\n", iq); }

        int32_t jqm = q->tdt - q->tbr; /* Max relative col of {j}'s trace. */
        for (int32_t j = 0; j < ncoq; j++) 
          { afq->uph[j] = (j <= jqm ? 0 : drtree_planar_NO_ROW); /* Measures down. */
            afq->loh[j] = (j <= jqm ? 0 : drtree_planar_NO_ROW); /* Measured up. */
          }

        /* Now place subtrees as specified by {abo,flp}: */
        for (int32_t r = 0; r < nch; r++)
          { int32_t ic = afq->chi[r];
            drtree_node_t *c = &(dt[ic]);
            drtree_planar_info_t *afc = &(af[ic]);
            afc->flp = flp[r];
            if (debug) { fprintf(stderr, "            r = %d ic = %d abo = %c flp = %c\n", r, ic, "FT"[abo[r]], "FT"[flp[r]]); }
            /* Stack subtree {c} above or below current tree of {q}: */
            int32_t rowr = 1; /* Min relative disp of {q}, unsigned. */ 
            int32_t ncoc = afc->jhi - afc->jlo + 1;
            for (int32_t jc = 0; jc < ncoc; jc++) 
              { /* Get proximal horizon of {c}'s tree: */
                int32_t ahcj = (abo[r] == flp[r] ? afc->uph[jc] : afc->loh[jc]);
                /* Get horizon of {q}'s tree towards {c}: */
                int32_t jq = jc + (c->tbr - q->tbr);
                assert((jq >= 0) && (jq < ncoq));
                int32_t bhqj = (abo[r] ? afq->uph[jq] : afq->loh[jq]);
                if (debug) { fprintf(stderr, "              jc = %d jq = %d ahcj = %d bhqj = %d\n", jc, jq, ahcj, bhqj); }
                /* Compute displacement of {r} rel to {q} to avoid overlap on col {jq}: */
                int32_t drj = ahcj + bhqj + 1;
                if (drj > rowr) { rowr = drj; }
              }
            /* Assign {c.row}: */
            afc->row = ( abo[r] ? +rowr : -rowr );
            if (debug) { fprintf(stderr, "            child row = %d\n", afc->row); }
            /* Update {q}'s horizons: */
            for (int32_t jc = 0; jc < ncoc; jc++) 
              { /* Get upper and lower horizons {uphcj,lohcj} of child's tree, WITH SIGNS: */
                int32_t uphcj = afc->row + (flp[r] ? afc->loh[jc] : afc->uph[jc]);
                int32_t lohcj = afc->row - (flp[r] ? afc->uph[jc] : afc->loh[jc]);
                if (debug) { fprintf(stderr, "              jc = %d child rows = {%d .. %d}", jc, uphcj, lohcj); }
                /* get upper and lower horizons {uphqj,lohqj} of {q}, WITH SIGNS: */
                int32_t jq = jc + (c->tbr - q->tbr);
                assert((jq >= 0) && (jq < ncoq));
                int32_t uphqj = + afq->uph[jq];
                int32_t lohqj = - afq->loh[jq];
                if (debug) { fprintf(stderr, " jq = %d rows before = {%d .. %d}", jq, uphqj, lohqj); }
                /* Check and update horizons of {q}: */
                if (abo[r])
                  { assert(lohcj > uphqj); }
                else
                  { assert(uphcj < lohqj); }
                if (uphcj > uphqj) { uphqj = uphcj; }
                if (lohcj < lohqj) { lohqj = lohcj; }
                /* Store new horizons, WITHOUT SIGNS: */
                afq->uph[jq] = + uphqj;
                afq->loh[jq] = - lohqj;
                if (debug) { fprintf(stderr, " after = {%d .. %d}\n", afq->uph[jq], afq->loh[jq]); }
              }
          }
        if (debug) { fprintf(stderr, "          < %s\n", __FUNCTION__); }
      }
      
     int32_t compute_envelope_radius(void)
       { 
         if (debug) { fprintf(stderr, "          > %s\n", __FUNCTION__); }
         if (debug) { fprintf(stderr, "            iq = %d\n", iq); }
         int32_t radius = 0;
         for (int32_t jq = 0; jq < ncoq; jq++)
           { assert((+ afq->uph[jq]) >= (- afq->loh[jq]));
             int32_t rhj = (int32_t)imax(abs(afq->uph[jq]), abs(afq->loh[jq]));
             if (rhj > radius) { radius = rhj; }
           }
         if (debug) { fprintf(stderr, "          < %s\n", __FUNCTION__); }
         return radius;
       }

    void free_child_horizons(void)
      { for (int32_t kc = 0; kc < nch; kc++)
          { int32_t ic = afq->chi[kc];
            assert((ic > iq) && (ic < ni));
            drtree_node_t *c = &(dt[ic]);
            drtree_planar_info_t *afc = &(af[ic]);
            assert(! drtree_node_is_null(c)); 
            assert(afc->uph != NULL); free(afc->uph); afc->uph = NULL;
            assert(afc->loh != NULL); free(afc->loh); afc->loh = NULL;
          }
      }
  }
  
int32_t drtree_planar_place_top_trees
  ( int32_t ni,
    drtree_node_t dt[],
    drtree_planar_info_t af[],
    int32_t nr,
    int32_t rut[],
    int32_t tMin,
    int32_t ncols,
    int32_t rdr[]
  )
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "      > %s\n", __FUNCTION__); }

    assert((ni >= 0) && (ni < drtree_indivs_MAX));

    /* Best order/flip arrangement found: */
    int32_t ord_best[nr]; /* Tells the order in which to try to stack the trees. */          
    bool_t flp_best[nr]; /* Tells which trees should be flipped up/down */         
    int32_t height_best = INT32_MAX; /* Envelope height for {abo_best,flp_best}. */ 
    
    int32_t uph[ncols]; /* Highest occ diagram row in each column. */
      
    auto void try_ord_flp(int32_t ord[], bool_t flp[]);
      /* To be given to {drtree_planar_enum_top_tree_placements}.
        Assumes {ord[0..nr-1]} is some permutation of the IDs of the
        root individuals, and {flp[0..nr-1]} are the corresponding
        up/down flip assignments. Computes the diagram height for that
        arrangement and updates {ord_best,flp_best,height_best}
        accordingly. */

    auto void assign_absolute_root_rows(int32_t ord[], bool_t flp[]);
      /* For each {r} in {0..nr-1}, sets the flip bit
        {af[iq].flp} of each root node {iq = ord[r]}. Then
        computes the absolute row displacement {af[iq].row} 
        by stacking the trees in the order specified by
        {ord[0..nr-1]}. Also computes the upper horizon {uph[0..ncols-1]} 
        of the diagram.
        
        Assumes that the horizons {af[iq].uph[0..q.nco-1]} and
        {af[iq].loh[0..q.nco-1]} of each root {q} are defined. */
    
    auto int32_t compute_diagram_height(void);
      /* Computes the height of the diagram, that is,
        the maximum occupied row plus one, based on the 
        diagram's upper horizon {uph[0..ncols-1]}. */

    drtree_planar_enum_top_tree_placements(nr, try_ord_flp);
     
    /* Set best placement, compute final horizons of {q}: */
    assign_absolute_root_rows(ord_best, flp_best);

    if (debug) { fprintf(stderr, "      < %s\n", __FUNCTION__); }
    return height_best;
    
    void try_ord_flp(int32_t ord[], bool_t flp[])
      { assign_absolute_root_rows(ord, flp);
        int32_t height = compute_diagram_height();
        if (height < height_best)
          { for (int32_t r = 0; r < nr; r++)
              { ord_best[r] = ord[r];
                flp_best[r] = flp[r];
                height_best = height; 
              }
          }
      }   
    
    void assign_absolute_root_rows(int32_t ord[], bool_t flp[])
      { 
        /* Initialize global horizon with {-1}: */
        for (int32_t j = 0; j < ncols; j++) { uph[j] = -1; }
        
        /* Now place top trees as specified by {ord,flp}: */
        for (int32_t r = 0; r < nr; r++)
          { int32_t iq = rut[ord[r]];
            drtree_node_t *q = &(dt[iq]);
            assert(! drtree_node_is_null(q));
            drtree_planar_info_t *afq = &(af[iq]);
            int32_t ncoq = afq->jhi - afq->jlo + 1;
            afq->flp = flp[r];
            /* Stack subtree {q} above previousy placed trees: */
            int32_t rowr = 0; /* Min absolute row of {q}. */ 
            for (int32_t jq = 0; jq < ncoq; jq++) 
              { /* Get lower horizon of {q}'s tree: */
                int32_t ahj = (flp[r] ? afq->uph[jq] : afq->loh[jq]);
                /* Get global upper horizon: */
                int32_t jg = jq + (q->tbr - tMin);
                assert((jg >= 0) && (jg < ncols));
                int32_t bhj = uph[jg];
                /* Compute row of {q} needed to avoid overlap on column {jg}: */
                int32_t drj = ahj + bhj + 1;
                if (drj > rowr) { rowr = drj; }
              }
            /* Assign {q.row}: */
            afq->row = rowr;
            /* Update the gobal horizon: */
            for (int32_t jq = 0; jq < ncoq; jq++) 
              { /* Get upper horizon of {q}'s tree: */
                int32_t ahj = afq->row + (flp[r] ? afq->loh[jq] : afq->uph[jq]);
                int32_t bhj = afq->row - (flp[r] ? afq->uph[jq] : afq->loh[jq]);
                /* Check and update global horizon: */
                int32_t jg = jq + (q->tbr - tMin);
                assert((jg >= 0) && (jg < ncols));
                assert(bhj > uph[jg]);
                if (ahj > uph[jg]) { uph[jg] = ahj; }
              }
          }
      }
      
    int32_t compute_diagram_height(void)
      { int32_t row_max = -1;
        for (int32_t j = 0; j < ncols; j++) 
          { if (uph[j] > row_max) { row_max = uph[j]; } }
        return row_max + 1;
      }
  }

void drtree_planar_set_absolute_rows
  ( int32_t ni,
    drtree_node_t dt[],
    drtree_planar_info_t af[],
    int32_t rdr[]
  )
  { 
    assert((ni >= 0) && (ni < drtree_indivs_MAX));
  
    /* Scan in chrono order (parent before children): */
    for (int32_t iq = 0; iq < ni; iq++)
      { drtree_node_t *q = &(dt[iq]);
        drtree_planar_info_t *afq = &(af[iq]);
        if (drtree_node_is_null(q))
          { rdr[iq] = -1; }
        else
          { int32_t ip = q->par;
            if (ip == -1)
              { /* Root -- {q.row} is absolute: */
                assert(afq->row >= 0);
                rdr[iq] = afq->row;
              }
            else
              { /* Not root -- {q.row} is relative to parent: */
                assert((ip >= 0) && (ip < iq));
                drtree_node_t *p = &(dt[ip]);
                drtree_planar_info_t *afp = &(af[ip]);
                assert(! drtree_node_is_null(p));
                if (afp->flp)
                  { rdr[iq] = rdr[ip] - afq->row; 
                    afq->flp = (! afq->flp);
                  }
                else
                  { rdr[iq] = rdr[ip] + afq->row; }
                assert(rdr[iq] >= 0);
              }
          }
      }
  }

void drtree_planar_info_free(int32_t ni, drtree_planar_info_t *af)
  { 
    assert((ni >= 0) && (ni < drtree_indivs_MAX));
    for (int32_t iq = 0; iq < ni; iq++)
      { drtree_planar_info_t *afq = &(af[iq]);
        if (afq->chi != NULL) { free(afq->chi); }
        if (afq->uph != NULL) { free(afq->uph); }
        if (afq->loh != NULL) { free(afq->loh); }
      }
    if (af != NULL) { free(af); }
  }

void drtree_planar_enum_subtree_placements
  ( int32_t nch,
    drtree_planar_abo_flp_proc_t try_abo_flp
  )
  { 
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "        > %s\n", __FUNCTION__); }
  
    if (debug) { fprintf(stderr, "          nch = %x\n", nch); }
            
    /* Tentative above/below and flip asignment: */
    bool_t abo[nch];
    bool_t flp[nch];

    /* Compute the number {nass} of possible {abo/flp} assignments, except
      for a global flip (reverse all {abo} and {flp}): */
    int32_t nass; /* Number of possible assignments, or INT32_MAX if too many: */
    if (nch > drtree_planar_enum_sub_num_MAX) 
      { nass = INT32_MAX; }
    else
      { nass = (nch == 0 ? 1 : (1 << (2*nch - 1))); }
    
    int32_t try_max = drtree_planar_enum_sub_try_MAX;
    if (nass > try_max)
      { /* Try only some arrangements: */
        /* Initialize with even children below, odd ones above: */
        for (int32_t r = 0; r < nch; r++) { abo[r] = ((r % 2) == 1); flp[r] = FALSE; }
        /* Iterate with random changes: */
        for (int32_t try = 0; try < try_max; try++)
          { try_abo_flp(abo, flp);
            /* Change a random {abo[ra]} or {flp[ra]}: */
            int32_t r = int32_abrandom(0,nch-1);
            if ((try % 2) == 0)
              { abo[r] = ! abo[r]; }
            else
              { flp[r] = ! flp[r]; }
          }
      } 
    else
      { /* Enumerate all combinations of {abo,flp} except for a global flip: */
        for (int32_t ass = 0; ass < nass; ass++)
          { /* Take {nch} bits of {ass} as {abo}, {nch-1} as {flp}: */
            if (debug) { fprintf(stderr, "            ass = %x\n", ass); }
            int32_t xx = ass;
            for (int32_t r = 0; r < nch; r++)
              { abo[r] = ((xx & 1) != 0);
                xx >>= 1;
                flp[r] = ((xx & 1) != 0);
                xx >>= 1;
              }
            if (nch > 0) { assert(! flp[nch-1]); }
            assert(xx == 0);
            try_abo_flp(abo, flp);
          }
      }
     if (debug) { fprintf(stderr, "        < %s\n", __FUNCTION__); }
  }

void drtree_planar_enum_top_tree_placements
  ( int32_t nr,
    drtree_planar_ord_flp_proc_t try_ord_flp
  )
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "        > %s\n", __FUNCTION__); }

    /* Tentative order and flip assignment: */
    int32_t ord[nr];
    bool_t flp[nr];

    /* Compute the number {nass} of possible order and flip assignemnts, except
      for a global flip (reverse all {ord} and {flp}): */
    int32_t nord; /* Number of possible orders, or INT32_MAX if too many: */
    int32_t nflp; /* Number of possible flips, or INT32_MAX if too many: */
    int32_t nass; /* Number of possible orders and flips, or INT32_MAX if too many: */
    if (nr > drtree_planar_enum_top_num_MAX) 
      { nord = INT32_MAX; nflp = INT32_MAX; nass = INT32_MAX; }
    else if (nr == 0)
      { nord = 1; nflp = 1; nass = 1; }
    else
      { nord = 1;
        for (int32_t k = 1; k <= nr; k++) { nord *= k; }
        nflp = (1 << (nr - 1));
        nass = nord*nflp;
      }
    if (debug) { fprintf(stderr, "          nr = %d nord = %d nflp = %d\n", nr, nord, nflp); }
      
    /* Initialize with identity perm, no flip: */
    for (int32_t r = 0; r < nr; r++) { ord[r] = r; flp[r] = FALSE; }
    
    int32_t try_max = drtree_planar_enum_top_try_MAX;
    if (nass > try_max)
      { /* Try only some arrangements, with random changes: */
        assert(nr >= 2);
        for (int32_t try = 0; try < try_max; try++)
          { if (debug)
              { fprintf(stderr, "            try = %d ord+flp = (", try);
                for (int32_t r = 0; r < nr; r++) { fprintf(stderr, " %d:%c", ord[r], "FT"[flp[r]]); }
                fprintf(stderr, " )\n");
              }
            try_ord_flp(ord, flp);
            /* Change a random {ord[ra]} or {flp[ra]}: */
            if ((try % 2) == 0)
              { int32_t r = int32_abrandom(0,nr-1);
                int32_t s = int32_abrandom(0,nr-2);
                if (s >= r) { s++; }
                int32_t t = ord[r]; ord[r] = ord[s]; ord[s] = t;
              }
            else
              { int32_t r = int32_abrandom(0,nr-1);
                flp[r] = ! flp[r];
              }
          }
      } 
    else
      { /* Enumerate all combinations of {ord,flp} except for a global flip: */
        for (int32_t ko = 0; ko < nord; ko++)
          { for (int32_t kf = 0; kf < nflp; kf++)
              { /* Take {nr-1} bits of {kf} as {flp}: */
                int32_t xx = kf;
                for (int32_t r = 0; r < nr-1; r++)
                  { flp[r] = ((xx & 1) != 0);
                    xx >>= 1;
                  }
                if (debug)
                  { fprintf(stderr, "            ko = %d kf = %d ord+flp = (", ko, kf);
                    for (int32_t r = 0; r < nr; r++) { fprintf(stderr, " %d:%c", ord[r], "FT"[flp[r]]); }
                    fprintf(stderr, " )\n");
                  }
                assert(! flp[nr-1]);
                try_ord_flp(ord, flp);
              }
            /* Find the decreasing tail: */
            int32_t kt = nr-1;
            while ((kt > 0) && (ord[kt-1] > ord[kt])) { kt--; }
            /* Reverse the decreasing tail: */
            int32_t j0 = kt, j1 = nr-1;
            while (j0 < j1) { int32_t t = ord[j0]; ord[j0] = ord[j1]; ord[j1] = t; j0++; j1--; }
            if (kt > 0)
              { /* Insert the next element {ord[kt-1]} in the order: */
                int32_t t = ord[kt-1];
                while ((kt < nr) && (ord[kt] < t)) { ord[kt-1] = ord[kt]; kt++; }
                ord[kt-1] = t;
              }
            else
              { assert(ko >= nord); }
          }
      }
    if (debug) { fprintf(stderr, "        < %s\n", __FUNCTION__); }
  }

void drtree_planar_check_info
  ( int32_t ni,
    drtree_node_t dt[],
    drtree_planar_info_t af[],
    int32_t tMin,
    int32_t tMax
  )
  {
    assert((ni >= 0) && (ni < drtree_indivs_MAX));
    assert(tMin < tMax);
    
    int32_t ncols = tMax - tMin + 1;
    
    /* Forward scan: check {jlo,jhi}: */
    for (int32_t iq = 0; iq < ni; iq++)
      { drtree_node_t *q = &(dt[iq]);
        drtree_planar_info_t *afq = &(af[iq]);
        if (drtree_node_is_null(q))
          { assert(afq->jlo > afq->jhi);
          }
        else
          { assert((afq->jlo >= 0) && (afq->jlo <= afq->jhi) && (afq->jhi < ncols));
            assert(afq->jlo == q->tbr - tMin);
            assert(afq->jhi >= q->tdt - tMin);
          }
        assert(afq->row == drtree_planar_NO_ROW);
        assert(afq->uph == NULL);        
        assert(afq->upMax == drtree_planar_NO_ROW); 
        assert(afq->loh == NULL);
        assert(afq->loMax == drtree_planar_NO_ROW);          
      }
      
    /* Reverse scan: check {nch,chi}: */
    int32_t nch[ni];
    for (int32_t iq = 0; iq < ni; iq++) { nch[iq] = 0; }
    for (int32_t iq = ni-1; iq >= 0; iq--)
      { drtree_node_t *q = &(dt[iq]);
        drtree_planar_info_t *afq = &(af[iq]);
        if (drtree_node_is_null(q))
          { assert(afq->nch == 0);
            assert(afq->chi == NULL);
          }
        else
          { assert(afq->nch == nch[iq]);
            assert((afq->nch == 0) == (afq->chi == NULL));
            /* Check children: */
            for (int32_t kc = 0; kc < afq->nch; kc++)
              { int32_t ic = afq->chi[kc];
                /* Chidlren must come after parent: */
                assert((ic > iq) && (ic < ni));
                drtree_node_t *c = &(dt[ic]);
                /* Chidlren can't be null: */
                assert(! drtree_node_is_null(c));
                /* Children must have {q} as parent: */
                assert(c->par == iq);
                /* Children must be in decreasing order of birth time: */
                if (kc > 0)
                  { int32_t ic1 = afq->chi[kc-1];
                    drtree_node_t *c1 = &(dt[ic1]);
                    assert(c1->tbr >= c->tbr);
                  }
              }
            /* Get parent, check life span, and bump its {nch}: */
            int32_t ip = q->par;
            if (ip != -1)
              { assert((ip >= 0) && (ip < iq));
                drtree_node_t *p = &(dt[ip]);
                drtree_planar_info_t *afp = &(af[ip]);
                /* Parent can't be null: */
                assert(! drtree_node_is_null(p));
                /* Parent's life span must include birth of {q}: */
                assert((afp->jlo <= afq->jlo) && (afp->jhi >=  afq->jlo));
                /* Bump its child count: */
                nch[ip]++;
              }
          }
      }
  }
