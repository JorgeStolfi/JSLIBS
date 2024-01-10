/* Last edited on 2011-12-24 02:28:06 by stolfilocal */

#include <delaunay.h>
#include <delaunay_plot.h>

#include <quad.h>
#include <bool.h>
#include <affirm.h>

#include <ioprotos.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

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
    
    double fnx = (double)nj; /* Number of steps in {x}. */
    double fny = (double)ni; /* Number of steps in {y}. */
    
    int i, j;
    int k = 0;
    if (polar) 
      { /* Polar grid with {ni} layers of {nj} sites, slightly twisted: */
        /* Angular increment: */
        double dt = 2*M_PI/(nj <= 1 ? 1.0 : nj-1);
        /* Radius reduction factor: */
        double rredf = 1.0 - (nj <= 12 ? 0.50 : (2*M_PI)/nj); 
        /* Relative twist between layers (fraction of a {j}-step): */
        double dtwist = 0.50/(ni <= 1 ? 1 : ni-1);
        double r = 1.0; /* Radius of current layer. */
        for (i=0; i < ni; i++) 
          { /* Add layer {i} at radius {r} counting from outside in: */
            for (j=0; j < nj; j++) 
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
      { /* Rectangular grid, slightly sheared and squeezed: */
         /* Compute grid steps: */
        double dy = 2.0/(ni <= 1 ? 1.0 : fny-1);
        double dx = 2.0/(nj <= 1 ? 1.0 : fnx-1);
        /* Relative X displacement between layers (fraction of a {j}-step): */
        double shear = 0.01;
        double shrink = 0.001;
        for (i=0; i < ni; i++) 
          { /* Add horizontal layer {i}: */
            double y = (i - fny/2)*dy;
            for (j=0; j < nj; j++) 
              { double x = ((j - fnx/2) + shear*(i - fny/2))*dx;
                double scale = 1.0 + shrink*(x*x + y*y); 
                st[k].index = k;
                r2_t p;
                p.c[0] = x/scale;
                p.c[1] = y/scale;
                st[k].pt = delaunay_hi2_from_r2(&p);
                k++;
              }
          }
      }
    if (VERBOSE) 
      { for (k = 0; k < nsites; k++)
          { assert(st[k].index == k);
            delaunay_debug_site ("given", &(st[k]));  fprintf(stderr, "\n");
            assert(abs(st[k].pt.c.c[0]) <= MC);
            assert(abs(st[k].pt.c.c[1]) <= MC);
            assert(abs(st[k].pt.c.c[2]) <= MC);
          }
      }
    (*nsitesp) = nsites; (*stp) = st;
  }
