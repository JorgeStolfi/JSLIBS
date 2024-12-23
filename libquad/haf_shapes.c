/* See haf_shapes.h. */
/* Last edited on 2024-12-22 10:31:04 by stolfi */

#define haf_shapes_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#include <stdint.h>
#include <stdlib.h>

#include <jslibs_copyright.h>
#include <bool.h>

#include <haf.h>
#include <haf_shapes.h>

haf_arc_t haf_shapes_torus(void)
  { haf_arc_t a, b;
    a = haf_make_stick(0); 
    b = haf_make_stick(1);
    haf_splice(a, b);
    haf_splice(haf_sym(a), a);
    haf_splice(haf_sym(b), a);
    return a;
  } 
     
haf_arc_t haf_shapes_star(uint32_t n)
  { haf_arc_t a, b;
    a = haf_make_stick(0);
    for(uint32_t k = 1; k < n; k++)
      { b = haf_make_stick(k);
        haf_splice(a, b);
        a = b;
      }
    return a;
  } 
     
haf_arc_t haf_shapes_ring(uint32_t n)
  { haf_arc_t fst, a, b;
    a = haf_make_stick(0); 
    fst = a;
    for (uint32_t k = 1;  k < n; k++)
      { b = haf_make_stick(k); 
        haf_splice(b, haf_sym(a));
        a = b;
      }
    haf_splice(fst, haf_sym(a));
    return fst;
  } 
     
haf_arc_t haf_shapes_orange(uint32_t n)
  { haf_arc_t a, b;
    a = haf_make_stick(0); 
    for (uint32_t k = 1;  k < n; k++)
      { b = haf_make_stick(k); 
        haf_splice(b, a);
        haf_splice(haf_sym(b), haf_sym(a));
        a = b;
      }
    return a;
  }

haf_arc_t haf_shapes_pyramid(uint32_t n)
  { /* Build an {n}-armed star: */
    haf_arc_t a = haf_shapes_star(n);
    /* Build an {n}-sided ring: */
    haf_arc_t b = haf_shapes_ring(n);
    /* Stitch them together: */
    haf_arc_t c = haf_sym(a);
    for(uint32_t k = 0; k < n; k++)
      { haf_splice(b, c);
        c = haf_dnext(c);
        b = haf_lnext(b);
      }
    return c;
  } 
  
haf_arc_t haf_shapes_tetra(void)
  { haf_arc_t s, t, a;
    s = haf_shapes_ring(4);
    a = haf_make_stick(5);
    t = haf_rnext(haf_rnext(s));
    
    haf_splice(s, a); 
    haf_splice(t, haf_sym(a));
    a = haf_make_stick(6);
    haf_splice(haf_sym(s), a);
    haf_splice(haf_sym(t), haf_sym(a));
    return a;
  }

haf_arc_t haf_shapes_cube (void)
  { haf_arc_t b, t, e;
    b = haf_shapes_ring(4);
    t = haf_shapes_ring(4); 
    t = haf_oprev(t);
    for (uint32_t k = 0;  k < 4; k++)
      { e = haf_make_stick(8+k);
        haf_splice(b, e); 
        haf_splice(haf_sym(e), t);
        b = haf_rprev(b); 
        t = haf_rnext(t);
      }
    return b;
  }

#define haf_shapes_C_AUTH \
  "J. Stolfi, IC-UNICAMP"
  
#define haf_shapes_C_HIST \
  "2023-10-05 J. Stolfi: Created."
