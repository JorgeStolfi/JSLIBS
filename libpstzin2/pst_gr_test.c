/* See {pst_gr_test.h} */
/* Last edited on 2025-03-15 14:17:37 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <jsprintf.h>
#include <jsrandom.h>

#include <pst_gr.h>
#include <pst_gr_path.h>

#include <pst_gr_test.h>

#define NONE pst_gr_NONE

pst_gr_t* pst_gr_test_make_graph(uint32_t nf, uint32_t nh)
  { bool_t debug = FALSE;
    demand(nf >= 1, "invalid {nf}");
    demand((nh >= 1) && ((nh % 2) == 1), "invalid {nh}");
    uint32_t nb = 7; /* Blocks. */
    uint32_t mx = nh + nf;     /* Additional cols per block. */
    uint32_t nx = nf + nb*mx;  /* Total cols. */
    uint32_t ny = 2*nf + nh;   /* Total rows.*/ 

    pst_gr_t *gr = pst_gr_new(200,200, nx,ny);
    
    if (debug) { fprintf(stderr, "  making the main grid...\n"); }
    pst_gr_vertex_t vtab[nx*ny]; /* Table of grid vertices, by row. */
    pst_gr_arc_t atab[nb];       /* One arc on side of each hole. */
    for (uint32_t y = 0; y < ny; y++)
      { for (uint32_t x = 0; x < nx; x++)
          { bool_t exists; /* True iff vertex {[x,y]} does exist. */
            int32_t bi;    /* Block index, or {-1}. */
            int32_t xrel;  /* Column relative to block, or {-1}. */
            if (x < nf)
              { xrel = -1; 
                bi = -1;
                exists = TRUE; 
              }
            else
              { xrel = (int32_t)((x - nf) % mx);
                bi = (int32_t)((x - nf) / mx);
                exists = ((y < nf) || (y >= ny-nf) || (xrel >= nh));  
              }
            pst_gr_vertex_t vxy = NONE;
            if (exists)
              { r2_t pxy = (r2_t){{ x, y }};
                vxy = pst_gr_add_vertex(gr, (int32_t)x, (int32_t)y, pxy);
                gr->vdata[vxy].vmark = (((x + y) % 2) == 0 ? 0 : 2);  
                if (debug) { fprintf(stderr, "  added vertex %d [%d,%d]\n", vxy, x, y); }
                for (int32_t edir = 0; edir <= 1; edir++)
                  { pst_gr_vertex_t dst;
                    if (edir == 0)
                      { dst = (x == 0 ? NONE : vtab[nx*y + x-1]); }
                    else
                      { dst = (y == 0 ? NONE : vtab[nx*(y-1) + x]); }
                    if (dst != NONE)
                      { double delta = dabrandom(-5,+5);
                        double weight = (drandom() < 0.10 ? 0.0 : dabrandom(0.10, 2.0));
                        pst_gr_path_t P = pst_gr_path_NULL;
                        if (debug) { fprintf(stderr, "  adding edge %d-%d\n", vxy, dst); }
                        pst_gr_arc_t axyd = pst_gr_add_edge(gr, vxy, dst, delta, weight, P, TRUE);
                        /* Remember one edge on the side of each hole: */
                        if ((bi >= 0) && (edir == 1) && (xrel == nh) && (y == nf))
                          { atab[bi] = pst_gr_arc_sym(axyd); }
                      }
                  }
              }
            vtab[nx*y + x] = vxy;
          }
      }
      
    if (debug) { fprintf(stderr, "  adding the special vertices...\n"); }
    for(uint32_t bi = 0; bi < nb; bi++)
      { uint32_t deg = bi; /* Degree of special vertex. */
        uint32_t xi = nf + bi*mx + nh/2;
        uint32_t yi = nf + nh/2;
        r2_t phi = (r2_t){{ xi, yi }};
        /* Add vertex in center of hole: */
        pst_gr_vertex_t vi = pst_gr_add_vertex(gr, (int32_t)xi, (int32_t)yi, phi);
        gr->vdata[vi].vmark = ((bi % 2) == 0 ? 1 : 4);  
        if (debug) { fprintf(stderr, "  added special vertex %d [%d,%d]\n", vi, xi, yi); }
        pst_gr_arc_t aj = atab[bi];
        uint32_t nah = 4*(nh + 1); /* Number of arcs in the hole. */
        uint32_t jah = 0;          /* Number of hole arcs already traversed. */
        for (int32_t j = 0; j < deg; j++)
          { assert(deg > 0);
            if (debug) { fprintf(stderr, "  aj = "); pst_gr_arc_print(stderr, gr, aj); }
            pst_gr_vertex_t vj = pst_gr_arc_org(gr, aj);
            gr->vdata[vj].vmark = 3;  
            r2_t *pj = &(gr->vdata[vj].coords);
            if (debug) { r2_gen_print(stderr, pj, "%8.4f", "  pj = ( ", " ", " )\n"); }
            double delta = dabrandom(-5,+5);
            double weight = (drandom() < 0.10 ? 0.0 : dabrandom(0.10, 2.0));
            pst_gr_path_t Pj = pst_gr_test_throw_path(pj, (uint32_t)j, &phi);
            if (debug) 
              { fprintf(stderr, "  adding edge %d-%d CCW of ", vi, vj);
                pst_gr_arc_print(stderr, gr, aj);
                r2_t ad = pst_gr_arc_start_dir(gr, aj);
                r2_gen_print(stderr, &ad, "%8.4f", " dir = ( ", " ", " )");
                fprintf(stderr, "\n"); 
              }
            pst_gr_arc_t rj = pst_gr_add_edge(gr, vj, vi, delta, weight, Pj, TRUE);
            if (debug) 
              { fprintf(stderr, "  added arc rj %d-%d = ", vi, vj); 
                pst_gr_arc_print(stderr, gr, rj);
                r2_t rd = pst_gr_arc_start_dir(gr, rj);
                r2_gen_print(stderr, &rd, "%8.4f", " dir = ( ", " ", " )");
                fprintf(stderr, "\n"); 
              }
            assert(pst_gr_arc_onext(gr, aj) == rj);
            /* Advance {aj} around the hole: */
            uint32_t kah = ((uint32_t)j + 1)*nah/deg;  /* Number of arcs in the hole. */
            while (jah < kah) { aj = pst_gr_arc_lnext(gr, aj); jah++; }
          }
      }
    return gr;
  }

pst_gr_path_t pst_gr_test_throw_path(r2_t *p0, uint32_t n, r2_t *p1)
  { pst_gr_path_t P = pst_gr_path_NULL;
    if (n > 0)
      { P.n = n;
        P.v = talloc(n, r2_t);
        /* Get a transversal vector {u}: */
        r2_t u; r2_sub(p1, p0, &u);
        r2_cross(&u, &u); 
        for (int32_t i = 0; i < n; i++)
          { double r = ((double)i+1)/(n+1);
            r2_t pi; r2_mix(1-r, p0, r, p1, &pi);
            double s = r*(1-r)*dabrandom(-0.5,+0.5);
            r2_mix(1.0, &pi, s, &u, &pi);
            P.v[i] = pi;
          }
      }
    return P;
  }
            
