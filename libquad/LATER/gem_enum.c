/* See gem_enum.h. */
/* Last edited on 2009-03-06 15:10:09 by stolfi */

#define gem_enum_C_copyright \
  "Copyright © 1996, 2006 Institute of Computing, Unicamp."

/* AUTHORS

  These enumeration procedures were originally developed as the
  Modula-3 file {Oct.m3} by J. Stolfi and Rober M. Rosi in 1993. The
  latter was converted to C by J. Stolfi in 1996 and was substantially
  revised by him in January 2007. */

#include <oct.h>
#include <gem_enum.h>
#include <malloc.h>
#include <stdint.h>

uint gem_node_degree(gem_arc_t p)
  { uint deg = 0;
    auto bool_t count_edge(gem_arc_t p);
    EnumIncidencesOfElem(Place_vec_desc(&p,1), 0, 1, &count_edge, NULL);
    return deg;
    bool_t count_edge(gem_arc_t p) { deg++; return FALSE; }
  }

uint DegreeOfWall(gem_arc_t p)
  { uint n = 0;
    gem_arc_t s = p;
    do { n++; s = NextF(s); } while (s != p);
    return n;
  }

uint DegreeOfEdge(gem_arc_t p)
  { uint n = 0;
    gem_arc_t s = p;
    do { n++; s = NextE(s); } while (s != p);
    return n;
  }

bool_t EnumPlaces(gem_arcs_t root, gem_visit_t visit, gem_arcs_t *vP)
  { 
    /* Pick three functions {stp[0..2]} that generate 
      all even walks. Then enumerate all places reachable
      by {<stp[0..2]>}. */
      
    gem_step_t stp[3];
    gem_visit_t vis[4];
    vis[0] = visit; stp[0] = &ONext; /* == Flip2.Flip1. */
    vis[1] = visit; stp[1] = &SymVF; /* == Flip0.Flip2. */
    vis[2] = visit; stp[2] = &NextF; /* == Flip3.Flip2. */
    vis[3] = visit;                   /* == next root. */
    
    return enum_orbits(root, 3, stp, vis, vP);
  }

/* ENUMERATING ELEMENTS */

bool_t EnumElemsEven(gem_arcs_t root, uint dimObj, gem_visit_t visit, gem_arcs_t *vP)
  { 
    demand(dimObj <= 3, "invalid dimension for target elements");
    /* Set {stp[0..2]} to flip pairs such that {stp[0..1]}
      generates the even subgroup of {< flp[i] : i != dimObj >}
      and {stp[2]} changes {p.elt[dimObj]}. So, each orbit of {stp[0..1]}
      corresponds to one evenly-reachable total orientation of one
      element of dimension {dimObj}, and vice-versa. */
    gem_step_t stp[3];
    switch(dimObj)
      {
        case 0: /* Nodes: */
          stp[0] = &Flip21; /* == Flip2.Flip1 = ONext. */
          stp[1] = &Flip32; /* == Flip3.Flip2 = NextF. */
          stp[2] = &Flip03; /* == Flip0.Flip3 = Clock. */
          break;

        case 1: /* Edges: */
          stp[0] = &Flip32; /* == Flip3.Flip2 = NextF. */
          stp[1] = &Flip03; /* == Flip0.Flip3 = Clock. */
          stp[2] = &Flip13; /* == Flip1.Flip3 = SymEC. */
          break;

        case 2: /* Walls: */
          stp[0] = &Flip10; /* == Flip1.Flip0 = PrevE. */
          stp[1] = &Flip03; /* == Flip0.Flip3 = Clock. */
          stp[2] = &Flip02; /* == Flip0.Flip2 = SymVF. */
          break;

        case 3: /* Cells: */
          stp[0] = &Flip21; /* == Flip2.Flip1 = ONext. */
          stp[1] = &Flip10; /* == Flip1.Flip0 = PrevE. */
          stp[2] = &Flip03; /* == Flip0.Flip3 = Clock. */
          break;

        default: assert(FALSE);
      }

    auto bool_t MyVisit(gem_arc_t p);
      /* Appends {p} to {vP} and calls {visitNode(p)}. */
    
    /* We should call {visit(p)} only after a {stp[3]} or a new root: */
    gem_visit_t vis[4];
    vis[0] = NULL;      
    vis[1] = NULL;      
    vis[2] = MyVisit;
    vis[3] = MyVisit;

    int nP = 0; /* Visited places are {vP->el[0..nP-1]}, if {vP != NULL}. */ 
    
    bool_t stop = enum_orbits(root, 3, stp, vis, NULL);
    
    /* Trim the visit list to the true size: */
    if (vP != NULL) { gem_arcs_trim(vP, nP); }
    
    /* Return the abort indicator: */
    return stop;
    
    bool_t MyVisit(gem_arc_t p)
      { if (vP != NULL) { Place_vec_expand(vP, nP); vP->el[nP] = p; nP++; }
        return (visit != NULL ? visit(p) : FALSE);
      }
  }

bool_t EnumElemsAny(gem_arcs_t root, uint dimObj, gem_visit_t visit, gem_arcs_t *vP)
  { demand(dimObj <= 3, "invalid dimension for target elements");
    /* Set {stp[0..3]} to a permutation of {flp[0..3]} so that
      {stp[3]==flp[dimObj]}. So, each orbit of {stp[0..1]} corresponds
      to one element of dimension {dimObj}, and vice-versa. */
    gem_step_t stp[4];
    stp[0] = &Flip0;
    stp[1] = &Flip1;
    stp[2] = &Flip2;
    stp[3] = &Flip3;
    /* Permute {stp[0..3]} so that {flp[dimObj]} is in {stp[3]}; */
    /* Try to preserve some sort of dual symmetry. */
    uint k = dimObj; /* Invariant: {stp[k] = flp[dimObj]}. */
    if (k < 2) 
      { /* Swap {stp[i]} with {stp[3-i]} for all {i}: */
        gem_step_t t;
        t = stp[0]; stp[0] = stp[3]; stp[3] = t;
        t = stp[1]; stp[1] = stp[2]; stp[2] = t;
        k = 3 - k;
      }
    /* Now the Invariant holds and {k} is {2} or {3}. */
    if (k == 2)
      { /* Swap even and odd elements of {stp[0..3]}: */
        gem_step_t t;
        t = stp[0]; stp[0] = stp[1]; stp[1] = t;
        t = stp[3]; stp[3] = stp[2]; stp[2] = t;
        k = 3;
      }
    /* Now the Invariant holds and {k} is {3}. */

    auto bool_t MyVisit(gem_arc_t p);
      /* Appends {p} to {vP} and calls {visitNode(p)}. */
    
    /* We should call {visit(p)} only after a {stp[3]} or a new root: */
    gem_visit_t vis[5];
    vis[0] = NULL;      
    vis[1] = NULL;      
    vis[2] = NULL;      
    vis[3] = MyVisit;
    vis[4] = MyVisit;

    int nP = 0; /* Visited places are {vP->el[0..nP-1]}, if {vP != NULL}. */ 
    
    bool_t stop = enum_orbits(root, 4, stp, vis, NULL);
    
    /* Trim the visit list to the true size: */
    if (vP != NULL) { gem_arcs_trim(vP, nP); }
    
    /* Return the abort indicator: */
    return stop;
    
    bool_t MyVisit(gem_arc_t p)
      { if (vP != NULL) { Place_vec_expand(vP, nP); vP->el[nP] = p; nP++; }
        return (visit != NULL ? visit(p) : FALSE);
      }
  }

/* ENUMERATING THE STAR AND BOUNDARY OF AN ELEMENT */

bool_t EnumPlacesOfElem(gem_arcs_t root, uint dimPiv, gem_visit_t visit, gem_arcs_t *vP)
  { 
    demand(dimPiv <= 3, "invalid dimension for pivot elements");

    /* We set {stp[0..1]} to generators of all the 
      even walks that do not change {p.elt[dim]}.
      Then we enumerate all places reachable through
      {<stp[0..1]>}. */
    gem_step_t stp[2];
    switch(dimPiv)
      {
        case 0: /* Node: */
          stp[0] = &Flip32; /* == Flip3.Flip2. */
          stp[1] = &Flip21; /* == Flip2.Flip1. */
          break;

        case 1: /* Edge */
          stp[0] = &Flip32; /* == Flip3.Flip2. */
          stp[1] = &Flip30; /* == Flip0.Flip3. */
          break;

        case 2: /* Wall */
          stp[0] = &Flip01; /* == Flip0.Flip1. */
          stp[1] = &Flip03; /* == Flip0.Flip3. */
          break;

        case 3: /* Cell */
          stp[0] = &Flip01; /* == Flip0.Flip1. */
          stp[1] = &Flip12; /* == Flip1.Flip2. */
          break;

        default: assert(FALSE);
      }
    
    auto bool_t MyVisit(gem_arc_t p);
      /* Appends {p} to {vP} and calls {visitNode(p)}. */
    
    /* We should visit {p} only after a {stp[2]} or a new root: */
    gem_visit_t vis[3];
    vis[0] = MyVisit;      
    vis[1] = MyVisit; 
    vis[2] = MyVisit; 
    
    int nP = 0; /* Visited places are {vP->el[0..nP-1]}, if {vP != NULL}. */ 
    
    /* OK, now enumerate the orbits: */
    bool_t stop = enum_orbits(root, 2, stp, vis, NULL);
    
    /* Trim the visit list to the true size: */
    if (vP != NULL) { gem_arcs_trim(vP, nP); }
    
    /* Return the abort indicator: */
    return stop;
    
    bool_t MyVisit(gem_arc_t p)
      { if (vP != NULL) { Place_vec_expand(vP, nP); vP->el[nP] = p; nP++; }
        return (visit != NULL ? visit(p) : FALSE);
      }
  }

bool_t EnumElemsOfElemEven(gem_arcs_t root, uint dimPiv, uint dimObj, gem_visit_t visit, gem_arcs_t *vP)
  { 
    demand(dimPiv <= 3, "invalid dimension for pivot elements");
    demand(dimObj <= 3, "invalid dimension for target elements");
    demand(dimPiv != dimObj, "pivot and target dimensions must be distinct");

    /* Set {stp[0..1]} to generators of all the even walks that do 
      not change {p.elt[dimPiv]}, so that only {stp[2]} changes
      {p.elt[dimObj]}. */
    gem_step_t stp[3];
    switch(dimPiv)
      {
        case 0: /* Pivot is node: */
          switch(dimObj)
            { case 1: /* edges */ stp[0] = &Flip23; stp[1] = &Flip31; break;
              case 2: /* walls */ stp[0] = &Flip31; stp[1] = &Flip12; break;
              case 3: /* cells */ stp[0] = &Flip12; stp[1] = &Flip23;break;
              default: assert(FALSE);
            }
          break;

        case 1: /* Pivot is edge: */
          switch(dimObj)
            { case 0: /* nodes */ stp[0] = &Flip32; stp[1] = &Flip30; break;
              case 2: /* walls */ stp[0] = &Flip03; stp[1] = &Flip02; break;
              case 3: /* cells */ stp[0] = &Flip20; stp[1] = &Flip23; break;
              default: assert(FALSE);
            }
          break;

        case 2: /* Pivot is wall: */
          switch(dimObj)
            { case 0: /* nodes */ stp[0] = &Flip31; stp[1] = &Flip01; break;
              case 1: /* edges */ stp[0] = &Flip03; stp[1] = &Flip13; break;
              case 3: /* cells */ stp[0] = &Flip10; stp[1] = &Flip30; break;
              default: assert(FALSE);
            }
          break;

        case 3: /* Pivot is cell: */
          switch(dimObj)
            { case 0: /* nodes */ stp[0] = &Flip12; stp[1] = &Flip01; break;
              case 1: /* edges */ stp[0] = &Flip20; stp[1] = &Flip12; break;
              case 2: /* walls */ stp[0] = &Flip01; stp[1] = &Flip20; break;
              default: assert(FALSE);
            }
          break;

        default: assert(FALSE);
      }

    /* We should visit {p} only after a {stp[2]} or a new root: */

    auto bool_t MyVisit(gem_arc_t p);
      /* Appends {p} to {vP} and calls {visitNode(p)}. */
    
    gem_visit_t vis[3];
    vis[0] = NULL;
    vis[1] = &MyVisit;
    vis[2] = &MyVisit;

    int nP = 0; /* Visited places are {vP->el[0..nP-1]}, if {vP != NULL}. */ 

    bool_t stop = enum_orbits(root, 2, stp, vis, vP);
    
    /* Trim the visit list to the true size: */
    if (vP != NULL) { gem_arcs_trim(vP, nP); }
    
    /* Return the abort indicator: */
    return stop;
    
    bool_t MyVisit(gem_arc_t p)
      { if (vP != NULL) { Place_vec_expand(vP, nP); vP->el[nP] = p; nP++; }
        return (visit != NULL ? visit(p) : FALSE);
      }
  }

bool_t EnumElemsOfElemAny(gem_arcs_t root, uint dimPiv, uint dimObj, gem_visit_t visit, gem_arcs_t *vP)
  { 
    demand(dimPiv <= 3, "invalid dimension for pivot elements");
    demand(dimObj <= 3, "invalid dimension for target elements");
    demand(dimPiv != dimObj, "pivot and target dimensions must be distinct");

    /* Set {stp[0..2]} to a permutation of {{ stp[i] : i != dimPiv }},
      with {flp[dimObj]} in {stp[2]}: */
    gem_step_t stp[4]; /* We will use {stp[3]} while sorting but not in the enumeration. */
    /* Permute {stp[0..3]} so that {flp[dimPiv]} is in {stp[3]}; */
    /* Try to preserve some sort of dual symmetry. */
    uint kPiv = dimPiv; /* Invariant 1: {stp[kPiv] = flp[dimPiv]}. */
    uint kObj = dimObj; /* Invariant 2: {stp[kObj] = flp[dimObj]}. */
    if (kPiv < 2) 
      { /* Swap {stp[i]} with {stp[3-i]} for all {i}: */
        gem_step_t t;
        t = stp[0]; stp[0] = stp[3]; stp[3] = t;
        t = stp[1]; stp[1] = stp[2]; stp[2] = t;
        kPiv = 3 - kPiv; kObj = 3 - kObj;
      }
    /* Now the Invariants 1,2 hold and {kPiv} is {2} or {3}. */
    if (kPiv == 2)
      { /* Swap even and odd elements of {stp[0..3]}: */
        gem_step_t t;
        t = stp[0]; stp[0] = stp[1]; stp[1] = t;
        t = stp[3]; stp[3] = stp[2]; stp[2] = t;
        kPiv = 3; kObj = (kObj ^ 1); 
      }
    /* Now the Invariants 1,2 hold and {kPiv} is {3}. */
    if (kObj == 0)
      { /* Swap {stp[0]} with {stp[2]} */
        gem_step_t t;
        t = stp[0]; stp[0] = stp[2]; stp[2] = t;
        kObj = 2; 
      }
    else if (kObj == 1)
      { /* Swap {stp[1]} with {stp[2]} */
        gem_step_t t;
        t = stp[1]; stp[1] = stp[2]; stp[2] = t;
        kObj = 2; 
      }
    /* Now the Invariants 1,2 hold and {kPiv} is {3} and {kObj} is 2. */

    /* We should visit {p} only after a {stp[2]} or a new root: */

    auto bool_t MyVisit(gem_arc_t p);
      /* Appends {p} to {vP} and calls {visitNode(p)}. */
    
    gem_visit_t vis[4];
    vis[0] = NULL;
    vis[1] = NULL;
    vis[2] = &MyVisit;
    vis[3] = &MyVisit;

    int nP = 0; /* Visited places are {vP->el[0..nP-1]}, if {vP != NULL}. */ 

    bool_t stop = enum_orbits(root, 2, stp, vis, vP);
    
    /* Trim the visit list to the true size: */
    if (vP != NULL) { gem_arcs_trim(vP, nP); }
    
    /* Return the abort indicator: */
    return stop;
    
    bool_t MyVisit(gem_arc_t p)
      { if (vP != NULL) { Place_vec_expand(vP, nP); vP->el[nP] = p; nP++; }
        return (visit != NULL ? visit(p) : FALSE);
      }
  }

bool_t EnumIncidencesOfElem(gem_arcs_t root, uint dimPiv, uint dimObj, gem_visit_t visit, gem_arcs_t *vP)
  { 
    demand(dimPiv <= 3, "invalid dimension for pivot elements");
    demand(dimObj <= 3, "invalid dimension for target elements");
    demand(dimPiv != dimObj, "pivot and target dimensions must be distinct");

    /* Get two {stp}s that geneerate all even walks that
      do not change {p.elt[dimPiv]}. The function {stp[0]}
      should also preserve {p.elt[dimObj]}. */
    gem_step_t stp[2];
    switch(dimPiv)
      {
        case 0: /* Pivot is node: */
          switch(dimObj)
            { case 1: /* edges */ stp[0] = &Flip23; stp[1] = &Flip13; break;
              case 2: /* walls */ stp[0] = &Flip13; stp[1] = &Flip23; break;
              case 3: /* cells */ stp[0] = &Flip12; stp[1] = &Flip32; break;
              default: assert(FALSE);
            }
          break;

        case 1: /* Pivot is edge: */
          switch(dimObj)
            { case 0: /* nodes */ stp[0] = &Flip23; stp[1] = &Flip03; break;
              case 2: /* walls */ stp[0] = &Flip30; stp[1] = &Flip32; break;
              case 3: /* cells */ stp[0] = &Flip02; stp[1] = &Flip32; break;
              default: assert(FALSE);
            }
          break;

        case 2: /* Pivot is wall: */
          switch(dimObj)
            { case 0: /* nodes */ stp[0] = &Flip13; stp[1] = &Flip10; break;
              case 1: /* edges */ stp[0] = &Flip30; stp[1] = &Flip10; break;
              case 3: /* cells */ stp[0] = &Flip01; stp[1] = &Flip03; break;
              default: assert(FALSE);
            }
          break;

        case 3: /* Pivot is cell: */
          switch(dimObj)
            { case 0: /* nodes */ stp[0] = &Flip12; stp[1] = &Flip10; break;
              case 1: /* edges */ stp[0] = &Flip02; stp[1] = &Flip01; break;
              case 2: /* walls */ stp[0] = &Flip01; stp[1] = &Flip02; break;
              default: assert(FALSE);
            }
          break;

        default: assert(FALSE);
      }

    /* We should visit {p} only after a {stp[2]} or a new root: */

    auto bool_t MyVisit(gem_arc_t p);
      /* Appends {p} to {vP} and calls {visitNode(p)}. */

    gem_visit_t vis[3];
    vis[0] = NULL;      
    vis[1] = MyVisit; 
    vis[2] = MyVisit; 

    int nP = 0; /* Visited places are {vP->el[0..nP-1]}, if {vP != NULL}. */ 

    /* OK, now enumerate the orbits: */
    bool_t stop = enum_orbits(root, 2, stp, vis, NULL);

    /* Trim the visit list to the true size: */
    if (vP != NULL) { gem_arcs_trim(vP, nP); }
    
    /* Return the abort indicator: */
    return stop;
    
    bool_t MyVisit(gem_arc_t p)
      { if (vP != NULL) { Place_vec_expand(vP, nP); vP->el[nP] = p; nP++; }
        return (visit != NULL ? visit(p) : FALSE);
      }
  }

/* ---------------------------------------------------------------------- */
/* JUNK? */

/* Enumerate edge quads */

void quad_do_enum (
    quad_arc_t a, 
    void visit_proc(quad_arc_t e, void *closure), 
    void *closure,
    unsigned mark
  );

unsigned next_mark = 1;

void quad_enum(
    quad_arc_t a, 
    void visit_proc(quad_arc_t e, void *closure), 
    void *closure
  )
  {
    unsigned mark = next_mark;
    next_mark++;
    if (next_mark == 0) next_mark = 1;
    quad_do_enum(a, visit_proc, closure, mark); 
  }

void quad_do_enum (
    quad_arc_t e, 
    void visit_proc(quad_arc_t e, void *closure), 
    void *closure,
    unsigned mark
  )
  {
    while (MARK(e) != mark)
      {
        visit_proc(e, closure);
        MARK(e) = mark;
        quad_do_enum (ONEXT(SYM(e)), visit_proc, closure, mark);
        e = ONEXT(e);
      }
  }

