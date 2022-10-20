/* Last edited on 2013-03-26 23:50:33 by stolfilocal */

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

void makesites(int ni, int nj, bool_t polar, int *nsitesp, delaunay_site_t **stp);

int main(int argc, char **argv);

int main(int argc, char **argv)
  { int ni = atoi(argv[1]);
    int nj = atoi(argv[2]);
    bool_t polar = (strcmp(argv[3], "polar") == 0);
    bool_t eps = (strcmp(argv[4], "eps") == 0);
    quad_arc_t e;
    delaunay_site_t *st = NULL;
    int nsites;
    fprintf(stderr, "Creating %d × %d sites...\n", ni, nj);
    makesites(ni, nj, polar, &nsites, &st);
    fprintf(stderr, "Building delaunay...\n");
    e = delaunay_build (st, nsites);
    fprintf(stderr, "Plotting delaunay...\n");
    plot_delaunay(e, st, nsites, "out/delgrid", eps);
    return(0);
  }

void makesites(int ni, int nj, bool_t polar, int *nsitesp, delaunay_site_t **stp)
  { 
    int nsites = ni*nj;
    delaunay_site_t *st = notnull(malloc(nsites*sizeof(delaunay_site_t)), "no mem");
    
    int i, j;
    int k = 0;
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
        for (i = 0; i < ni; i++) 
          { /* Add layer {i} at radius {r} counting from outside in: */
            for (j = 0; j < nj; j++) 
              { double t = (j + i*dtwist)*dt;
                st[k].index = k;
                r2_t p;
                p.c[0] = r*cos(t);
                p.c[1] = r*sin(t);
                st[k].pt = delaunay_hi2_from_r2(&p);
                k++;
              }
            r = rredf * r;
          }
      }
    else 
      { /* Rectangular grid, slightly sheared: */
        /* Compute grid parameters: */
        int32_t wx = 1, xmin = 0, dx = 0; /* {xmin,dx} divided by {wx}. */
        int32_t wy = 1, ymin = 0, dy = 0; /* {ymin,dy} divided by {wy}. */
        int32_t dsx = 0;                  /* Shear step {dsx} divided by {wx}. */
        if (nj > 1) 
          { wx = nj-1; xmin = -(nj-1); dx = 2; }
        if (ni > 1) 
          { wy = ni-1; ymin = -(ni-1); dy = 2;
            wx = wx*nj; xmin = xmin*nj; dx = dx*nj; dsx = 1;
          }
        for (i = 0; i < ni; i++) 
          { /* Add horizontal layer {i}: */
            int32_t y = ymin + i*dy; /* Divided by {wy} */
            for (j = 0; j < nj; j++) 
              { double x = xmin + j*dx + i*dsx; /* Divided by {wx} */
                st[k].index = k;
                hi2_point_t pt;
                pt.c.c[0] = wx*wy;
                pt.c.c[1] = wy*x;
                pt.c.c[2] = wx*y;
                st[k].pt = pt;
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
