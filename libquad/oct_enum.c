/* See oct_enum.h. */
/* Last edited on 2025-01-09 23:21:53 by stolfi */

#define oct_enum_C_copyright \
  "Copyright © 1996, 2006 State University of Campinas (UNICAMP).\n\n" jslibs_copyright ""

/* These enumeration procedures were originally developed as the
  Modula-3 file {Oct.m3} by J. Stolfi and Rober M. Rosi in 1993. The
  latter was converted to C by J. Stolfi in 1996 and was substantially
  revised by him in January 2007. */

#include <malloc.h>
#include <stdint.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <enum_orbits.h>

#include <oct.h>
#include <oct_enum.h>

bool_t oct_enum
  ( oct_arc_vec_t root, 
    uint NS, 
    oct_step_t *step[], 
    oct_visit_t *visit[], 
    oct_arc_vec_t *vP
  )
  { assert(sizeof(ref_t) == sizeof(oct_arc_t));
    ref_vec_t *root_x = (ref_vec_t *)&(root);
    enum_step_t **step_x = (enum_step_t **)step;
    enum_visit_t **visit_x = (enum_visit_t **)visit;
    ref_vec_t *vP_x = (ref_vec_t *)vP;
    return enum_items(*root_x, NS, step_x, visit_x, vP_x);
  }

bool_t oct_enum_cycle
  ( oct_arc_t root, 
    oct_step_t *step, 
    oct_visit_t *visit, 
    oct_arc_vec_t *vP
  )
  {
    assert(sizeof(ref_t) == sizeof(oct_arc_t));
    ref_t root_x = (ref_t)root;
    enum_step_t *step_x = (enum_step_t *)step;
    enum_visit_t *visit_x = (enum_visit_t *)visit;
    ref_vec_t *vP_x = (ref_vec_t *)vP;
    return enum_cycle(root_x, step_x, visit_x, vP_x);
  }

bool_t oct_enum_orbits
  ( oct_arc_vec_t root,
    uint NI, 
    oct_step_t *istep[], 
    uint NO, 
    oct_step_t *ostep[], 
    oct_visit_t *visit,
    oct_arc_vec_t *vP
  )
  {
    uint32_t nvis; /* Counts arcs added to {vP}. */
    
    auto bool_t do_visit(oct_arc_t p);
    /* A visit-func to be called when starting a new root or 
       a new {istep}-orbit.  It appends {p} to {vP} at position
       {nvis}, calls the client-given visit-func, and increments
       {nvis}. */
    
    /* Create the step-funcs and visit-funcs for enumeration: */
    uint32_t NS = NI + NO;
    oct_step_t *stp[NS];
    oct_visit_t *vis[NS+1];
    for (uint32_t i = 0; i < NS; i++)
      { stp[i] = (i < NI ? istep[i] : ostep[i - NI]);
        vis[i] = (i < NI ? NULL : &do_visit);
      }
    vis[NS] = &do_visit;
    
    nvis = 0;
    bool_t res = oct_enum(root, NS, stp, vis, NULL);
    if (vP != NULL) { oct_arc_vec_trim(vP, nvis); }
    return res;
    
    /* IMPEMENTATION OF LOCAL FUNCTIONS */
    
    bool_t do_visit(oct_arc_t p)
      { if (vP != NULL) 
          { oct_arc_vec_expand(vP, (int32_t)nvis);
            vP->e[nvis] = p;
          }
        nvis++;
        return (visit == NULL ? FALSE : visit(p));
      }
  }

bool_t oct_enum_arcs(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Step functions. */
    uint32_t NS = 2;
    oct_step_t *stp[NS];
    stp[0] = &oct_sym;
    stp[2] = &oct_onext;
    assert(NS == 2);
    
    /* Visit functions. */
    oct_visit_t *vis[NS+1];
    vis[0] = visit;
    vis[1] = visit;
    vis[2] = visit;
    assert(NS+1 == 3);
    
    return oct_enum(root, NS, stp, vis, vP);
  }

bool_t oct_enum_nodes(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Inner step functions: cycle around origin node. */
    uint32_t NI = 1;
    oct_step_t *istep[NI];
    istep[0] = &oct_onext;
    assert(NI == 2);
    
    /* Outer step functions: turn arc 180 degrees. */
    uint32_t NO = 1;
    oct_step_t *ostep[NO];
    ostep[0] = &oct_sym;
    assert(NO == 1);
    
    return oct_enum_orbits(root, NI, istep, NO, ostep, visit, vP);
  }

bool_t oct_enum_edges(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Inner step functions: turn arc 180 degrees. */
    uint32_t NI = 1;
    oct_step_t *istep[NI];
    istep[0] = &oct_sym;
    assert(NI == 1);
    
    /* Outer step functions: cycle around origin node. */
    uint32_t NO = 1;
    oct_step_t *ostep[NO];
    ostep[0] = &oct_onext;
    assert(NO == 1);
    
    return oct_enum_orbits(root, NI, istep, NO, ostep, visit, vP);
  }

bool_t oct_enum_faces(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Inner step functions: cycle around left face. */
    uint32_t NI = 1;
    oct_step_t *istep[NI];
    istep[0] = &oct_lnext;
    assert(NI == 1);
    
    /* Outer step functions: turn arc around. */
    uint32_t NO = 1;
    oct_step_t *ostep[NO];
    ostep[0] = &oct_sym;
    assert(NO == 1);
    
    return oct_enum_orbits(root, NI, istep, NO, ostep, visit, vP);
  }

bool_t oct_enum_octets(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Inner step functions: all tumblings. */
    uint32_t NI = 2;
    oct_step_t *istep[NI];
    istep[0] = &oct_rot;
    istep[1] = &oct_fflip;
    assert(NI == 2);
    
    /* Outer step functions: cycle around node. */
    uint32_t NO = 1;
    oct_step_t *ostep[NO];
    ostep[0] = &oct_onext;
    assert(NO == 1);
    
    return oct_enum_orbits(root, NI, istep, NO, ostep, visit, vP);
  }
