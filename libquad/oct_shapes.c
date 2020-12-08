/* See oct_shapes.h. */
/* Last edited on 2011-12-22 15:02:49 by stolfilocal */

#define oct_shapes_C_copyright \
  "Copyright © 1996, 2006 Institute of Computing, Unicamp."

#include <oct_shapes.h>

#include <oct.h>
#include <bool.h>

#define _GNU_SOURCE
#include <stdlib.h>

oct_arc_t make_ring(int n)
  { oct_arc_t fst, a, b;
    a = oct_make_edge(); 
    fst = a;
    uint i;
    for (i = 1; i < n; i++)
      { b = oct_make_edge(); 
        oct_splice(b, oct_sym(a));
        a = b;
      }
    oct_splice(fst, oct_sym(a));
    return fst;
  }

oct_arc_t make_orange(uint n)
  { return oct_tor(make_ring(n)); }

oct_arc_t make_tetra(void)
  { oct_arc_t s, t, a;
    s = make_ring(4);
    a = oct_make_edge();
    t = oct_rnext(oct_rnext(s));
    
    oct_splice(s, a); 
    oct_splice(t, oct_sym(a));
    a = oct_make_edge();
    oct_splice(oct_sym(s), a);
    oct_splice(oct_sym(t), oct_sym(a));
    return a;
  }

oct_arc_t make_stick(void)
  { oct_arc_t a;
    a = oct_make_edge();
    oct_splice(a, oct_sym(a));
    return a;
  }

oct_arc_t make_cube (void)
  { oct_arc_t b, t, e;
    b = make_ring(4);
    t = make_ring(4); 
    t = oct_oprev(t);
    uint i;
    for (i = 0; i < 4; i++)
      { e = oct_make_edge();
        oct_splice(b, e); 
        oct_splice(oct_sym(e), t);
        b = oct_rprev(b); 
        t = oct_rnext(t);
      }
    return b;
  }

oct_arc_t make_sausage (uint len)
  { oct_arc_t t;
    t = make_ring(2);
    buld_tower(2, len, t);
    return t;
  }

oct_arc_t make_fork(uint np, uint len)
  { oct_arc_t o;
    o = make_orange(np);
    uint i;
    for (i = 0; i < np; i++)
      { o = oct_onext(o); 
        buld_tower(2, len, o);
      }
    return o;
  }
  
void buld_tower(uint m, uint h, oct_arc_t a)
  { oct_arc_t s, e, t;
    t = a;
    uint i;
    for (i = 0; i < h; i++)
      { t = oct_oprev(t);
        s = make_ring(m);
        uint j;
        for (j = 0; j < m; j++)
          { e = oct_make_edge();
            oct_splice(t, e); oct_splice(oct_sym(e), s);
            t = oct_lnext(t); s = oct_rnext(s);
          }
        t = s;
      }
  }
  
oct_arc_t make_projetive (void)
  { oct_arc_t e;
    e = oct_make_edge();
    oct_splice (oct_fflip(oct_sym(e)), e);
    return e;
  }

oct_arc_t make_torus(void)
  { oct_arc_t a, b;
    a = oct_make_edge();
    b = oct_make_edge();
    oct_splice(a, b);
    oct_splice(oct_sym(a), a);
    oct_splice(oct_sym(b), a);
    return a;
  }
 
oct_arc_t make_bitorus(void)
  { oct_arc_t a, b, c, d, e, f;
    a = oct_make_edge();
    b = oct_make_edge();
    c = oct_make_edge();
    d = oct_make_edge();
    e = oct_make_edge();
    f = oct_make_edge();
    oct_splice(a, b);
    oct_splice(b, c);
    oct_splice(c, d);
    oct_splice(d, e);
    oct_splice(e, oct_sym(c));
    oct_splice(oct_sym(a), oct_sym(e));
    oct_splice(oct_sym(e), oct_sym(f));
    oct_splice(oct_sym(f), oct_sym(d));
    oct_splice(oct_sym(d), oct_sym(b)); 
    oct_splice(oct_sym(b), f); 
    return a;
  }
 
oct_arc_t make_tritorus(void)
  {
    oct_arc_t a[8];
    oct_arc_t t, s, e;
    /* Build two tetrahedra, inner and outer. Store in {a[i, k]} one
      arc such that {oct_left(a[i*4+k])} is face {k} of tetrahedron {i}. */
    uint i;
    for (i = 0; i <= 1; i++)
      { a[i*4+0] = make_tetra();
        a[i*4+1] = oct_onext(a[i*4+0]);
        a[i*4+2] = oct_sym(oct_onext(oct_sym(a[i*4+0])));
        a[i*4+3] = oct_oprev(a[i*4+2]);
      }
    /* Build tubes connecting corresponding faces of the two tetrahedra.
      Each tube has three edges and three square faces. */
    uint k;
    for (k = 0; k <= 3; k++)
      { t = a[0*4+k]; s = a[1*4+k];
        uint j;
        for (j = 1; j <= 3; j++)
          { e = oct_make_edge();
            oct_splice(t, e); 
            oct_splice(s, oct_fflip(oct_sym(e)));
            t = oct_lnext(t); s = oct_lnext(s);
          }
      }
    return s;
  }
 
oct_arc_t make_klein(void)
  { oct_arc_t a, b;
    a = oct_make_edge(); 
    b = oct_make_edge();
    oct_splice(a, b);
    oct_splice(oct_sym(a), a);
    oct_splice(oct_fflip(oct_sym(b)), a);
    return a;
  }
  
oct_arc_t make_klein2(void)
  { oct_arc_t a, b, c, d;
    
    a = oct_make_edge(); 
    b = oct_make_edge();
    c = oct_make_edge(); 
    d = oct_make_edge();
    
    /* Vertex 1: */
    oct_splice(b, a);
    oct_splice(a, oct_fflip(oct_sym(c)));
    oct_splice(oct_fflip(oct_sym(c)), oct_sym(a));
    
    /* Vertex 2: */
    oct_splice(c, d);
    oct_splice(d, oct_sym(b));
    oct_splice(oct_sym(b), oct_sym(d));

    return a;
  }
  
oct_arc_t make_klein3(void)
  { oct_arc_t a, b, c, d, e, f;
    a = oct_make_edge(); 
    b = oct_make_edge();
    c = oct_make_edge(); 
    d = oct_make_edge();
    e = oct_make_edge();
    f = oct_make_edge();
    
    /* Vertex 1: */
    oct_splice(b, a);
    oct_splice(a, oct_fflip(oct_sym(c)));
    oct_splice(oct_fflip(oct_sym(c)), oct_sym(a));
    
    /* Vertex 2: */
    oct_splice(f, d);
    oct_splice(d, oct_sym(b));
    oct_splice(oct_sym(b), oct_sym(d));
    
    /* Vertex 3: */
    oct_splice(c, e);
    oct_splice(e, oct_sym(f));
    oct_splice(oct_sym(f), oct_sym(e));
    
    return a;
  }
  
#define oct_shapes_C_AUTH \
  "Modula-3 version {MakeShape.m3} created by Rober M. Rosi"\
  " and J. Stolfi, ca. 1994."
  
#define oct_shapes_C_HIST \
  "Converted to C by J. Stolfi in fev/2007."
