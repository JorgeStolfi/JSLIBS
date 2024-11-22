/* Last edited on 2023-02-20 06:22:48 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <quad.h>
#include <bool.h>
#include <affirm.h>
#include <r2.h>
#include <hi2.h>

#include <delaunay.h>
#include <delaunay_plot.h>
#include <delaunay_debug.h>


#define VERBOSE 0

#define MC (delaunay_MAX_COORD)

void makesites(int32_t ni, int32_t nj, bool_t polar, int32_t *nsitesp, delaunay_site_t **stp);

int32_t main(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
  { int32_t ni = atoi(argv[1]);
    int32_t nj = atoi(argv[2]);
    char *gridType = argv[3];
    bool_t polar = (strcmp(gridType, "polar") == 0);
    
    delaunay_site_t *st = NULL;
    int32_t nsites;
    fprintf(stderr, "Creating %d × %d sites...\n", ni, nj);
    makesites(ni, nj, polar, &nsites, &st);

    fprintf(stderr, "Building delaunay...\n");
    quad_arc_t e = delaunay_build (st, nsites);

    fprintf(stderr, "Plotting delaunay...\n");
    
    char *prefix = "out/delgrid";
    char *tag = gridType;
    int32_t capLines = 2;
    epswr_figure_t *eps = delaunay_plot_new_figure(prefix, tag, capLines, st, nsites); 
    
    epswr_text(eps, "Delaunay and Voronoi", FALSE, 0.5, TRUE, FALSE);
    char *capt2 = jsprintf("%s grid sites", gridType);
    epswr_text(eps, capt2, FALSE, 0.5, TRUE, FALSE);
    free(capt2);
        
    delaunay_plot_diagram(eps, e, st, nsites);
    epswr_end_figure(eps);
    return(0);
  }

void makesites(int32_t ni, int32_t nj, bool_t polar, int32_t *nsitesp, delaunay_site_t **stp)
  { 
    int32_t nsites = ni*nj;
    delaunay_site_t *st = notnull(malloc(nsites*sizeof(delaunay_site_t)), "no mem");
    
    int32_t k = 0;
    if (polar) 
      { /* Polar grid with {ni} layers of {nj} sites, slightly twisted: */
        /* Angular increment: */
        double dt = 2*M_PI/nj;
        /* Radius of first and last ring: */
        double rmax = 1.0;
        double rmin = (ni <= 1 ? rmax : 0.05);
        /* Radius reduction factor from ring to next: */
        double rredf = (ni <= 1 ? 1.0 : pow(rmin/rmax, 1.0/(ni - 1))); 
        /* Relative twist between layers (fraction of a {j}-step): */
        double dtwist = 1.0/ni;
        double r = rmax; /* Radius of current layer. */
        for (int32_t i = 0; i < ni; i++) 
          { /* Add layer {i} at radius {r} counting from outside in: */
            for (int32_t j = 0; j < nj; j++) 
              { double t = (j + i*dtwist)*dt;
                delaunay_site_t *stk = &(st[k]);
                stk->index = k;
                r2_t p;
                p.c[0] = r*cos(t);
                p.c[1] = r*sin(t);
                stk->pt = delaunay_hi2_from_r2(&p);
                k++;
              }
            r = rredf * r;
          }
      }
    else 
      { /* Rectangular grid, slightly sheared: */
        /* Compute grid parameters: */
        int64_t wx, xmin, dx; /* {xmin,dx} divided by {wx}. */
        int64_t wy, ymin, dy; /* {ymin,dy} divided by {wy}. */
        int64_t dsx;      /* Shear step {dsx} divided by {wx}. */
        if (nj > 1) 
          { wx = nj-1; xmin = -(nj-1); dx = 2; }
        else 
          { wx = 1; xmin = 0; dx = 0; }
        if (ni > 1) 
          { wy = ni-1; ymin = -(ni-1); dy = 2;
            wx = wx*nj; xmin = xmin*nj; dx = dx*nj; dsx = 1;
          }
        else
          { wy = 1; ymin = 0; dy = 0; dsx = 0; }
        int64_t c[3];
        for (int32_t i = 0; i < ni; i++) 
          { /* Add horizontal layer {i}: */
            int64_t y = ymin + i*dy; /* Divided by {wy} */
            for (int32_t j = 0; j < nj; j++) 
              { int64_t x = xmin + j*dx + i*dsx; /* Divided by {wx} */
                delaunay_site_t *stk = &(st[k]);
                stk->index = k;
                c[0] = wx*wy;
                c[1] = wy*x;
                c[2] = wx*y;
                for (int32_t ic = 0; ic < 3; ic++)
                  { affirm(c[ic] < +(int64_t)MC, "overflow in coordinate");
                    affirm(c[ic] > -(int64_t)MC, "underflow in coordinate");
                    stk->pt.c.c[ic] = (int32_t)c[ic];
                  }
                k++;
              }
          }
      }
    if (VERBOSE) 
      { for (k = 0; k < nsites; k++)
          { assert(st[k].index == k);
            deldebug_print_site ("given", &(st[k]));  fprintf(stderr, "\n");
            assert(abs(st[k].pt.c.c[0]) <= MC);
            assert(abs(st[k].pt.c.c[1]) <= MC);
            assert(abs(st[k].pt.c.c[2]) <= MC);
          }
      }
    (*nsitesp) = nsites; (*stp) = st;
  }
