/* See {psgr_test.h} */
/* Last edited on 2025-04-27 03:23:05 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <jsprintf.h>
#include <jsrandom.h>
#include <r2.h>
#include <i2.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_path.h>

#include <psgr_test.h>

/* SHORTER NAMES */  

#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE

#define CLEAR   psgr_mark_CLEAR
#define KILL    psgr_mark_KILL   
#define KEEP    psgr_mark_KEEP   
#define DELETED psgr_mark_DELETED

psgr_t* psgr_test_make_graph(uint32_t nf, uint32_t nh)
  { bool_t debug = FALSE;
    demand(nf >= 1, "invalid {nf}");
    demand((nh >= 1) && ((nh % 2) == 1), "invalid {nh}");
    
    /* Each hole "eats" {nh} by {nh} vertices. */
    /* A block consists of a hole and frames at top, bootom, and right. */
    /* A vertical frame has {nf} by {2*nf + nh} vertices. */
    /* A horizontal frame has {nh} by {nf} vertices. */
    uint32_t nb = 7; /* Blocks. */
    uint32_t mx = nh + nf;     /* Additional cols per block. */
    uint32_t nx = nf + nb*mx;  /* Total cols. */
    uint32_t ny = 2*nf + nh;   /* Total rows.*/ 

    psgr_t *gr = psgr_new(200,200, nx,ny);
    
    if (debug) { fprintf(stderr, "  making the main grid...\n"); }
    psgr_vertex_t vtab[nx*ny]; /* Table of grid vertices, by row. */
    psgr_arc_t atab[nb];       /* One arc on side of each hole. */
    for (uint32_t y = 0; y < ny; y++)
      { for (uint32_t x = 0; x < nx; x++)
          { bool_t exists; /* True iff vertex {[x,y]} does exist. */
            int32_t bi;    /* Block index, or {-1}. */
            int32_t xrel;  /* Column relative to block, or {-1}. */
            if (x < nf)
              { /* In the leftmost vertical frame: */
                xrel = -1; 
                bi = -1;
                exists = TRUE; 
              }
            else
              { /* Inside some block: */
                xrel = (int32_t)((x - nf) % mx);
                bi = (int32_t)((x - nf) / mx);
                exists = ((y < nf) || (y >= ny-nf) || (xrel >= nh));  
              }
            psgr_vertex_t vxy = VNONE;
            if (exists)
              { i2_t xy = (i2_t){{ (int32_t)x, (int32_t)y }};
                vxy = psgr_add_vertex(gr, &xy);
                psgr_mark_t mkxy = CLEAR;
                psgr_vertex_set_mark(gr, vxy, mkxy);  
                if (debug) { fprintf(stderr, "  added vertex %d [%d,%d]\n", vxy.ix, x, y); }
                for (int32_t edir = 0; edir <= 1; edir++)
                  { psgr_vertex_t dst;
                    if (edir == 0)
                      { dst = (x == 0 ? VNONE : vtab[nx*y + x-1]); }
                    else
                      { dst = (y == 0 ? VNONE : vtab[nx*(y-1) + x]); }
                    if (dst.ix != VNONE.ix)
                      { double delta = dabrandom(-5,+5);
                        double weight = dabrandom(0.01, 2.0);
                        psgr_path_t P = psgr_path_NULL;
                        if (debug) { fprintf(stderr, "  adding edge v%d--v%d\n", vxy.ix, dst.ix); }
                        psgr_arc_t axyd = psgr_add_edge(gr, vxy, dst, delta, weight, P);
                        /* Remember one edge on the side of each hole: */
                        if ((bi >= 0) && (edir == 1) && (xrel == nh) && (y == nf))
                          { atab[bi] = psgr_arc_sym(axyd); }
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
        /* Add vertex in center of hole: */
        i2_t ipi = (i2_t){{ (int32_t)xi, (int32_t)yi }};
        psgr_vertex_t vi = psgr_add_vertex(gr, &ipi);
        psgr_mark_t mki = (psgr_mark_t[]){ CLEAR, KILL, KEEP }[bi % 3];
        assert(mki != DELETED);
        psgr_vertex_set_mark(gr, vi, mki); 
        if (debug) { fprintf(stderr, "  added special vertex %d [%d,%d]\n", vi.ix, xi, yi); }
        psgr_arc_t aj = atab[bi];
        uint32_t nah = 4*(nh + 1); /* Number of arcs in the hole. */
        uint32_t jah = 0;          /* Number of hole arcs already traversed. */
        for (int32_t j = 0; j < deg; j++)
          { assert(deg > 0);
            if (debug) { fprintf(stderr, "  aj = "); psgr_arc_print(stderr, gr, aj); }
            psgr_vertex_t vj = psgr_arc_org(gr, aj);
            psgr_mark_t mkj = CLEAR;
            if (mki == DELETED)
              { mkj = KILL; }
            else if (mki == CLEAR)
              { mkj = KEEP; }
            psgr_vertex_set_mark(gr, vj, mkj);  
            i2_t ipj = psgr_vertex_pixel(gr, vj);
            if (debug) { i2_gen_print(stderr, &ipj, "%d", "  ipj = [ ", " ", " ]\n"); }
            double delta = dabrandom(-5,+5);
            double weight = (drandom() < 0.10 ? 0.0 : dabrandom(0.10, 2.0));
            psgr_path_t Pj = psgr_path_throw(&ipi, (uint32_t)j, &ipj);
            if (debug) 
              { fprintf(stderr, "  adding edge v%d-v%d to origin of ", vi.ix, vj.ix);
                psgr_arc_print(stderr, gr, aj);
                r2_t ad = psgr_arc_starting_dir(gr, aj);
                r2_gen_print(stderr, &ad, "%8.4f", " dir = ( ", " ", " )");
                fprintf(stderr, "\n"); 
              }
            psgr_arc_t rj = psgr_add_edge(gr, vi, vj, delta, weight, Pj);
            psgr_edge_set_mark(gr, psgr_arc_edge(rj), CLEAR);
            if (debug) 
              { fprintf(stderr, "  added arc rj v%d-v%d = ", vi.ix, vj.ix); 
                psgr_arc_print(stderr, gr, rj);
                r2_t rd = psgr_arc_starting_dir(gr, rj);
                r2_gen_print(stderr, &rd, "%8.4f", " dir = ( ", " ", " )");
                fprintf(stderr, "\n"); 
              }
            assert(psgr_arc_onext(gr, aj).ix == psgr_arc_sym(rj).ix);
            /* Advance {aj} around the hole: */
            uint32_t kah = ((uint32_t)j + 1)*nah/deg;  /* Number of arcs in the hole. */
            while (jah < kah) { aj = psgr_arc_lnext(gr, aj); jah++; }
          }
      }
    return gr;
  }
            
