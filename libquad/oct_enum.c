/* See oct_enum.h. */
/* Last edited on 2023-10-05 12:18:06 by stolfi */

#define oct_enum_C_copyright \
  "Copyright © 1996, 2006 State University of Campinas (UNICAMP).\n\n" jslibs_copyright"

/* These enumeration procedures were originally developed as the
  Modula-3 file {Oct.m3} by J. Stolfi and Rober M. Rosi in 1993. The
  latter was converted to C by J. Stolfi in 1996 and was substantially
  revised by him in January 2007. */

#define _GNU_SOURCE
#include <malloc.h>
#include <stdint.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <enum_orbits.h>

#include <oct.h>
#include <oct_enum.h>

bool_t oct_enum
  ( oct_arc_vec_t root, 
    uint ns, 
    oct_step_t *step[], 
    oct_visit_t *visit[], 
    oct_arc_vec_t *vP
  )
  { assert(sizeof(ref_t) == sizeof(oct_arc_t));
    ref_vec_t *root_x = (ref_vec_t *)&(root);
    enum_step_t **step_x = (enum_step_t **)step;
    enum_visit_t **visit_x = (enum_visit_t **)visit;
    ref_vec_t *vP_x = (ref_vec_t *)vP;
    return enum_items(*root_x, ns, step_x, visit_x, vP_x);
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
    uint ni, 
    oct_step_t *istep[], 
    uint no, 
    oct_step_t *ostep[], 
    oct_visit_t *visit,
    oct_arc_vec_t *vP
  )
  {
    int32_t nvis; /* Counts arcs added to {vP}. */
    
    auto bool_t do_visit(oct_arc_t p);
    /* A visit-func to be called when starting a new root or 
       a new {istep}-orbit.  It appends {p} to {vP} at position
       {nvis}, calls the client-given visit-func, and increments
       {nvis}. */
    
    /* Create the step-funcs and visit-funcs for enumeration: */
    int32_t ns = ni + no;
    oct_step_t *stp[ns];
    oct_visit_t *vis[ns+1];
    int32_t i;
    for (i = 0; i < ns; i++)
      { stp[i] = (i < ni ? istep[i] : ostep[i - ni]);
        vis[i] = (i < ni ? NULL : &do_visit);
      }
    vis[ns] = &do_visit;
    
    nvis = 0;
    bool_t res = oct_enum(root, ns, stp, vis, NULL);
    if (vP != NULL) { oct_arc_vec_trim(vP, nvis); }
    return res;
    
    /* IMPEMENTATION OF LOCAL FUNCTIONS */
    
    bool_t do_visit(oct_arc_t p)
      { if (vP != NULL) 
          { oct_arc_vec_expand(vP, nvis);
            vP->e[nvis] = p;
          }
        nvis++;
        return (visit == NULL ? FALSE : visit(p));
      }
  }

bool_t oct_enum_arcs(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Step functions. */
    int32_t ns = 2;
    oct_step_t *stp[ns];
    stp[0] = &oct_sym;
    stp[2] = &oct_onext;
    assert(ns == 2);
    
    /* Visit functions. */
    oct_visit_t *vis[ns+1];
    vis[0] = visit;
    vis[1] = visit;
    vis[2] = visit;
    assert(ns+1 == 3);
    
    return oct_enum(root, ns, stp, vis, vP);
  }

bool_t oct_enum_nodes(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Inner step functions: cycle around origin node. */
    int32_t ni = 1;
    oct_step_t *istep[ni];
    istep[0] = &oct_onext;
    assert(ni == 2);
    
    /* Outer step functions: turn arc 180 degrees. */
    int32_t no = 1;
    oct_step_t *ostep[no];
    ostep[0] = &oct_sym;
    assert(no == 1);
    
    return oct_enum_orbits(root, ni, istep, no, ostep, visit, vP);
  }

bool_t oct_enum_edges(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Inner step functions: turn arc 180 degrees. */
    int32_t ni = 1;
    oct_step_t *istep[ni];
    istep[0] = &oct_sym;
    assert(ni == 1);
    
    /* Outer step functions: cycle around origin node. */
    int32_t no = 1;
    oct_step_t *ostep[no];
    ostep[0] = &oct_onext;
    assert(no == 1);
    
    return oct_enum_orbits(root, ni, istep, no, ostep, visit, vP);
  }

bool_t oct_enum_faces(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Inner step functions: cycle around left face. */
    int32_t ni = 1;
    oct_step_t *istep[ni];
    istep[0] = &oct_lnext;
    assert(ni == 1);
    
    /* Outer step functions: turn arc around. */
    int32_t no = 1;
    oct_step_t *ostep[no];
    ostep[0] = &oct_sym;
    assert(no == 1);
    
    return oct_enum_orbits(root, ni, istep, no, ostep, visit, vP);
  }

bool_t oct_enum_octets(oct_arc_vec_t root, oct_visit_t *visit, oct_arc_vec_t *vP)
  {  
    /* Inner step functions: all tumblings. */
    int32_t ni = 2;
    oct_step_t *istep[ni];
    istep[0] = &oct_rot;
    istep[1] = &oct_fflip;
    assert(ni == 2);
    
    /* Outer step functions: cycle around node. */
    int32_t no = 1;
    oct_step_t *ostep[no];
    ostep[0] = &oct_onext;
    assert(no == 1);
    
    return oct_enum_orbits(root, ni, istep, no, ostep, visit, vP);
  }
