#ifndef intgprod_H
#define intgprod_H

/* Integrals-of-products */
/* Last edited on 2009-03-05 10:41:27 by stolfi */

#define _GNU_SOURCE
#include <interval.h>

/* 
  GOAL
  
  This interface provides tools for computing the definite integral of
  a function {F(x)} that is the product of {n} functions
  {f[0..n-1](x)} over a given interval {X} of {x}-values. The
  functions {f[i]} must be non-negative; otherwise they may be
  arbitrary.
  
  However, the algorithm is optimized for situations where most of the
  {f[i]} have a few tall and narrow peaks separated by fairy flat
  areas. An example is the likelyhood integral that arises when
  fitting a one-parameter model to {n} observations, where the error
  distributions have non-negligible long-tailed distributions.
  
  GRAPH ENCLOSURES
  
  The procedure uses the concept of a /graph enclosure/ for a real
  function {f}, which is a set of pairs {(x,y)} that contains the 
  points {(x,f(x))} for all {x} within some region {X} of the domain.
  
  The graph enclosure used here is the Cartesian product of the domain
  subset {X} and an interval {[yLO _ yHI]} of values, such that {yLO
  <= f(x) <= yHI} for all {x} in {X}. If {X} itself is an interval,
  then the enclosure is a /box/, an axis-aligned rectangle of the X-Y
  plane, that contains the graph of {f} within that X-interval. Such
  box enclosures can be easily computed by interval arithmetic (IA). */

/*
  FACTOR PRESENTATION

  The factors {f[i]} are actually given by a procedure {est_f(i,X,Y,xm)}
  which takes a box enclosure {X × Y} of {f[i](x)} over a given
  interval {X} of {x} values, and an abscissa {xm} in {X}; and returns
  two intervals {Y1,Y2} that enclose {f[i](x)} over the two
  sub-intervals {X1=[LO(X),xm]} and {X2=[xm,HI(X)]}. */
  
typedef void intgprod_est_func_t
  ( int i,          /* Index of desired factor. */
    interval_t *X,  /* An interval of X values. */
    interval_t *Y,  /* Interval enclosure of {f[i](x)} for {x} in {*X}. */
    double xm,      /* The split point (interior to {*X}). */
    interval_t *Y1, /* (OUT) Enclosure for {f[i](x)} in {*X1}. */
    interval_t *Y2  /* (OUT) Enclosure for {f[i](x)} in {*X2}. */
  );

/* 
  BRANCH-AND-BOUND INTEGRATION
  
  The following algorithm computes the integral of the product {F(x)}
  by a branch-and-bound technique. The state of the computation is a
  partition of the integration range {X} into some number {m} of
  intervals {X[0..m-1]} (the /X-cells/), and interval enclosures for
  the values of {F(x)} and/or the factors {f[i](x)} when {x} ranges
  over each of these X intervals. The boundaries between the X-cells
  (the /X-vertices/) are here denoted {x[0..m]} where {x[0] = LO(X)}
  and {x[m] = HI(X)}.
  
  This information is stored in a binary tree {T}. Each node {t} of
  {T} is associated to an X-interval {t.X}, which is the union of one
  or more consecutive X-cells {X[ilo..ihi]}; and contains an interval
  {t.IoF}, that is guaranteed to enclose the definite integral of
  {F(x)} over the interval {t.X}.
  
  In addition, a node {t} may be either a /leaf node/ or a /split
  node/. A split node contains a /splitting abscissa/ {t.xm}, and
  pointers to two (non-null) children nodes {t0,t1}, which have
  {HI(t0.X) == xm == LO(t1.X)}. A leaf node (which has {t.X} reduced to a
  single X-cell {X[k]}) contains instead information {t.nv,t.iv,t.fY,t,fP}
  on the ranges of the factors {f[i](x)} for {x} in that X-cell.
  
  Note that the interval {t.X} does not have to be stored, since it
  can be computed on the fly as one traverses the tree {T}. */
 
typedef struct intgprod_node_t 
  { interval_t IoF; /* Enclosure for the integral of {F} over {t.X}. */
    
    /* For split nodes: */
    double xm;                /* Splitting X-coordinate, or {NAN} if leaf node. */
    struct intgprod_node_t *ch[2]; /* Children of {t}, or {NULL} if leaf node. */
    
    /* For leaf nodes: */
    int nv;         /* Number of factors that are bounded separately. */
    int *iv;        /* The indices of the variable factors are {t.iv[0..nv-1]}. */
    interval_t *fY; /* The intervals {t.Y[0..nv-1]} are enclosures for those factors. */
    interval_t fP;  /* Enclosure for the product of the remaining factors. */
  } intgprod_node_t;
  /* Type of a node {t} of the tree that describes the state of the computation.  */
   
intgprod_node_t *intgprod_node_new(void);
  /* Allocates a new {intgprod_node_t} record {t}, with {t.IoF = empty}.
    {t.ch={NULL,NULL}}, {t.xm = NAN}, and {t.Y = NULL}. */

void intgprod_node_free(intgprod_node_t *t);
  /* Recaims the area used by {t}, including its {t.Y} vector
    if not NULL. */
   
/* 
  FACTOR RANGE DATA
  
  For each leaf node {t}, there is a set {t.V} of {nv} /variable/
  factors, whose indices are {t.iv[0..nv-1]}. For these factors, the
  algorithm keeps separate interval enclosures {t.fY[0..nv-1]}. The
  algorithm also keeps a single interval {t.fP} that encloses the
  product of all other factors not in this set.
  
  The goal is to keep in {t.fd.V} only those factors that may have
  significant variation withing the interval {t.X}. Keeping separate
  enclosures for those factors makes it easier to recompute the
  enclosures when the leaf has to be split. Note that if {f[i]} is
  known to be monotonic in the X coordinate, then its enclosures
  {Y1,Y2} in the children intervals {X1,X2} of {t.X} of can often be
  computed just from the enclosure {t.Y[j]}, where {iv[j] = i}, and
  the value {ym = f[i](xm)}.
  
  On the other hand, if a factor has little relative variation within
  the range {t.X}, then one does not need to call {est_f} again when
  splitting the leaf --- since the same Y interval can be used for
  both halves. In that case, the enclosures of all those factors can
  be multiplied together. */

void intgprod_condense_factor_list(intgprod_node_t *t, double minVar);
  /* Scans the set {t.fY[0..t.nv-1]} of separate factor bounds of a leaf node {t},
    and removes any entry whose relative variation (the width of its interval enclosure
    divided by its lower bound) is less than {minVar}.  These intervals
    are multiplied into the product interval {t.fP}. Also updates
    the index list {t.iv[0..t.nv-1]} and the count {t.nv}. */
  
/*
  BRANCH-AND-BOUND INTEGRATION
  
  The desired integral is obviously contained in the interval {r.IoF}
  of the root node {r}. If that interval is not sufficiently precise,
  we can refine it by splitting some leaf node {t}. That means
  introducing a new split abscissa {xm} in {t.X} that breaks it into
  two sub-cells {X0} and {X1}; and calling {est_f} to obtain the
  factor enclosure information for the two new children {t0,t1} of {t}.
  
  Note also that the interval {t.IoF} of a split node is simply the
  sum of the {.IoF} intervals of its two children, in the sense of
  interval arithmetic (IA); whereas the {t.IoF} field of a leaf node
  is simply the width of {t.X} times the interval enclosure of the
  product of the enclosures {t.Z[0..n-1]}. */  

void intgprod_node_split
  ( intgprod_node_t *t, 
    interval_t *X, 
    intgprod_est_func_t *est_f
  );
  /* Assumes that {*X} is the X-span of {t}.  If {t} is a leaf node,
    makes it into a split node by choosing a point {xm} somewhere
    inside {*X}.  If {t} is already a split node, calls itself recursively on
    the child that contributes the most to the undertainty of {t.IoF}.
    In either case, updates the fields of {t} to preserve the
    invariants. */

interval_t intgprod_integrate(double a[]);
  /* Assumes that all factors are monotonic in the interval. */
  /* We know F[i](a[0]) and F[i](a[1]) for all {i}. */
  
  /* Pick a value {amd} in {(a[0] _ a[1])} */
  /* Evaluate {F[i](amd)} for all {i}. */
  /* Get an interval for the integral in each half. */
  /* !!! FINISH !!! */

#endif
