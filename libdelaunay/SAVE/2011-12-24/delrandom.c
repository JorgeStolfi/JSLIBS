/* Last edited on 2011-12-24 01:05:02 by stolfilocal */

#include <delaunay.h>
#include <delaunay_plot.h>

#include <quad.h>
#include <affirm.h>
#include <pswr.h>
#include <jsrandom.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define NORMAL 1
#define NTOSS 9
#define VERBOSE 1
#define RBIG 10.0

#define MC (delaunay_MAX_COORD)

delaunay_site_t *makesites(int nsites, bool_t normal);

int main(int argc, char **argv);

int main(int argc, char **argv)
  { int nsites = atoi(argv[1]);
    bool_t normal = (strcmp(argv[2], "normal") == 0);
    bool_t eps = (strcmp(argv[3], "eps") == 0);
    quad_arc_t e;
    fprintf(stderr, "Creating %d sites...\n", nsites);
    delaunay_site_t *st = makesites(nsites, normal);
    fprintf(stderr, "Building delaunay...\n");
    e = delaunay_build (st, nsites);
    fprintf(stderr, "Plotting delaunay...\n");
    plot_delaunay(e, st, nsites, "out/delrandom", eps);
    return(0);
  }

delaunay_site_t *makesites(int nsites, bool_t normal)
  { delaunay_site_t *st = notnull(malloc(nsites*sizeof(delaunay_site_t)), "no mem");
    
    srandom(4615);

    int k;
    for (k=0; k < nsites; k++) 
      { st[k].index = k;
        /* Generate a random point {p} in the square {[-1 _ =1]^2}: */
        r2_t p;
        if (normal) 
          { /* Central limit approx to normal distribuition: */
            double s = 0, t = 0;
            int j;
            for (j=0; j<NTOSS; j++)
              { s += 2 * (drandom() - 0.5);
                t += 2 * (drandom() - 0.5);
              }
            p.c[0] = s/sqrt(NTOSS);
            p.c[1] = t/sqrt(NTOSS);
          }
        else 
          { /* Uniform distribution: */
            p.c[0] = 2 * (drandom() - 0.5);
            p.c[1] = 2 * (drandom() - 0.5);
          }
        /* Convert {p} to integer homogeneous coordinates: */
        st[k].pt = delaunay_hi2_from_r2(&p);
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
    return st;
  }
