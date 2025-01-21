/* See {haf.h}. */
/* Last edited on 2025-01-10 08:14:11 by stolfi */
 
#define haf_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright
 
/* Written by J. Stolfi in October 2023, loosely parallel to {quad.c}
  but without the {rot} (dual) operators. */ 

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <affirm.h>

#include <haf.h>

typedef struct half_rep_t { } half_rep_t; 
  /* Fictitious target of a {haf_arc_t} pointer. */

typedef struct haf_edge_rec_t 
  { haf_arc_t lnext[2]; /* Next edge counterclockwise around left face. */
    haf_edge_id_t eid;                  /* Sequential edge number. */
  } haf_edge_rec_t;

vec_typeimpl(haf_arc_vec_t, haf_arc_vec, haf_arc_t);
vec_typeimpl(haf_edge_vec_t, haf_edge_vec, haf_edge_t);
vec_typeimpl(haf_edge_id_vec_t, haf_edge_id_vec, haf_edge_id_t);
vec_typeimpl(haf_arc_id_vec_t, haf_arc_id_vec, haf_arc_id_t);

/* BIT HACKING OF ARC ADDRESSES */

#define EDGE(a) ((haf_edge_t)(((uint64_t)(a)) & (~(uint64_t)1)))
  /* Pointer to the {haf_edge_t} record underlying the arc pointer {a}.
    CANNOT be assigned to. */

#define EDGE_ID(a) ((EDGE(a))->eid)
  /* ID number of the edge underlying the arc pointer {a}. CAN be assigned to. */

#define DIR_BIT(a) (((uint64_t)(a)) & ((uint64_t)1))
  /* The direction bit of arc pointer {a}: 0 for the base arc, 
    1 for its opposite. CANNOT be assigned to. */

#define ARC_ID(a) (2*EDGE_ID(a) + DIR_BIT(a))
  /* ID number of the arc pointer {a}. CANNOT be assigned to. */

#define SYM(a) ((haf_arc_t )(((uint64_t)(a)) ^ ((uint64_t)1)))
  /* The {a.sym} arc pointer operator through bit hacking. CANNOT be assigned to. */

#define LNEXT(a) ((EDGE(a))->lnext[DIR_BIT(a)])
  /* The {a.lnext} arc pointer operator through bit hacking. CAN be assigned to. */

#define ARC(e,db) ((haf_arc_t )(((uint64_t)(e)) + ((uint64_t)(db))))
  /* The arc pointer with underlying edge record pointer {e} and direction bit {db}.
    CANNOT be assigned to. */

/* IMPLEMENTATIONS */

haf_arc_t haf_sym(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return SYM(a);
  }
  
haf_arc_t haf_lnext(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return LNEXT(a); }
 
haf_arc_t haf_rnext(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return SYM(LNEXT(SYM(a))); }
       
haf_arc_t haf_oprev(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return LNEXT(SYM(a)); }
  
haf_arc_t haf_dprev(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return SYM(LNEXT(a)); }

haf_arc_t haf_lprev(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    haf_arc_t c = a;
    while (TRUE)
      { haf_arc_t b = LNEXT(c);
        if (b == a) { break; }
        c = b;
      }
    return c;
  }
     
haf_arc_t haf_dnext(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    haf_arc_t c = a;
    while (TRUE)
      { haf_arc_t b = SYM(LNEXT(c));
        if (b == a) {break; }
        c = b;
      }
    return c;
  }

haf_arc_t haf_rprev(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return SYM(haf_lprev(SYM(a)));
  }
    
haf_arc_t haf_onext(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return SYM(haf_dnext(SYM(a)));
  }

haf_edge_t haf_edge(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return EDGE(a);
  }

haf_dir_bit_t haf_dir_bit(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return DIR_BIT(a);
  }

haf_arc_t haf_orient(haf_edge_t e, haf_dir_bit_t db)
  { demand(e != NULL, "invalid {NULL} argument");
    return ARC(e,db);
  }

haf_edge_id_t haf_edge_id(haf_edge_t e)
  { demand(e != NULL, "invalid {NULL} argument");
    return EDGE_ID(e);
  }

haf_arc_id_t haf_arc_id(haf_arc_t a)
  { demand(a != NULL, "invalid {NULL} argument");
    return ARC_ID(a);
  }

haf_arc_t haf_base_arc(haf_edge_t e)
  { demand(e != NULL, "invalid {NULL} argument");
    return ARC(e,0);
  }

haf_arc_t haf_make_stick(haf_edge_id_t eid)
  {
    demand(eid <= haf_edge_id_MAX, "edge identifier too big");
    
    haf_edge_rec_t *e = talloc(1, haf_edge_rec_t);
    e->eid = eid;
    
    haf_arc_t a = ARC(e,0);
    haf_arc_t b = ARC(e,1);
    
    e->lnext[0] = b;
    e->lnext[1] = a;

    return a;
  }

haf_arc_t haf_make_loop(haf_edge_id_t eid)
  {
    demand(eid <= haf_edge_id_MAX, "edge identifier too big");
    
    haf_edge_rec_t *e = talloc(1, haf_edge_rec_t);
    e->eid = eid;
    
    haf_arc_t a = ARC(e,0);
    haf_arc_t b = ARC(e,1);
    
    e->lnext[0] = a;
    e->lnext[1] = b;

    return a;
  }
  
void haf_edge_free(haf_edge_t e)
  { haf_arc_t a = ARC(e,0);
    haf_arc_t b = SYM(a);
    
    demand((e->lnext[0] == a) || (e->lnext[0] == b), "edge is not isolated");
    demand((e->lnext[1] == a) || (e->lnext[1] == b), "edge is not isolated");

    free(e);
  }

void haf_splice(haf_arc_t a, haf_arc_t b)
  { demand(a != NULL, "invalid {NULL} argument {a}");
    demand(b != NULL, "invalid {NULL} argument {b}");
    
    haf_arc_t an = LNEXT(a);
    haf_arc_t bn = LNEXT(b);
    LNEXT(a) = bn;
    LNEXT(b) = an;
  }

void haf_set_lnext(haf_arc_t a, haf_arc_t b)
  { demand(a != NULL, "invalid {NULL} argument {a}");
    demand(b != NULL, "invalid {NULL} argument {b}");
    LNEXT(a) = b;
  }

void haf_set_edge_id(haf_edge_t e, haf_edge_id_t eid)
  { 
    demand(e != NULL, "invalid {NULL} argument {a}");
    demand(eid <= haf_edge_id_MAX, "edge identifier too big");
    EDGE_ID(e) = eid;
  }

void haf_check_topology(haf_edge_count_t NE, haf_edge_t e[], haf_edge_id_t eid0, bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> %s\n", __FUNCTION__); } 
        
    demand(NE <= haf_edge_count_MAX, "too many edges");
    demand(eid0 <= haf_edge_id_MAX - NE + 1, "edge identifiers too big");
    
    haf_arc_count_t NA = 2*NE;   /* Number of arcs in structure and size of {cid}. */
    haf_arc_id_t aid0 = 2*eid0;  /* Lowest arc ID. */

    /* Check that the edge ids match indices in {a}: */
    if (verbose) { fprintf(stderr, "  checking edge IDs...\n"); } 
    for (haf_edge_count_t ke = 0; ke < NE; ke++) 
      { haf_edge_t edk = e[ke];
        demand(edk != NULL, "given arc is {NULL}");
        haf_edge_id_t edk_id = EDGE_ID(edk);
        demand(edk_id == ke + eid0, "inconsistent arc and edge ids");
        if (verbose) { fprintf(stderr, "    haf_edge_id(e[%lu]) = %lu\n", ke, edk_id); }
        haf_arc_t ak = ARC(edk,0);
        haf_arc_id_t ak_id = ARC_ID(ak);
        demand(ak_id == 2*edk_id, "inconsistent arc and edge ids");
        haf_dir_bit_t ak_dir = DIR_BIT(ak);
        if (verbose) { fprintf(stderr, "    haf_base_arc(e[%lu]) = arc[%lu] = %lu:%d\n", ke, 2*ke, edk_id, ak_dir); }
        demand(edk_id == eid0 + ke, "edge id does not match index");
      }
    
    /* No need to check that {.sym} is involution since the data structure ensures it. */
   
    /* Check that the {.lnext} pointers are arc permutations: */
    if (verbose) { fprintf(stderr, "  checking {.lnext} pointers...\n"); } 
    haf_arc_t *lprev = talloc(NA, haf_arc_t); /* Inverse of {.lnext}. */
    for (haf_arc_count_t ia = 0; ia < NA; ia++) { lprev[ia] = NULL; }

    auto void check_lnext(haf_arc_t c);
      /* Uses and updates the {lprev} vector. */

    for (haf_edge_count_t ke = 0; ke < NE; ke++) 
      { haf_edge_t edk = e[ke];
        haf_arc_t ak = ARC(edk,0);
        check_lnext(ak);
        check_lnext(SYM(ak));
      }
      
    if (verbose) { fprintf(stderr, "< %s\n", __FUNCTION__); } 
    return;
      
    void check_lnext(haf_arc_t b)
      { assert(b != NULL); /* Already checked, but just in case... */
        haf_arc_t u = LNEXT(b);
        haf_edge_id_t ueid = EDGE_ID(b);
        demand(ueid <= haf_edge_id_MAX, "invalid edge id");
        demand((ueid >= eid0) && (ueid < eid0 + NE), "edge id out of range");
        haf_arc_count_t ui = ARC_ID(u) - aid0;
        assert(ui < NA);
        demand(lprev[ui] == NULL, "{.lnext} is not a permutation");
        lprev[ui] = b;
      }
  }
