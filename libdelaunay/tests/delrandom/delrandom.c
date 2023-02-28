/* Last edited on 2023-02-20 06:22:57 by stolfi */

#include <delaunay.h>
#include <delaunay_plot.h>
#include <delaunay_debug.h>

#include <quad.h>
#include <affirm.h>
#include <epswr.h>
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

delaunay_site_t *makesites(int32_t nsites, bool_t normal);

int32_t main(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
  { int32_t nsites = atoi(argv[1]);
    char *distr = argv[2];
    bool_t normal = (strcmp(distr, "normal") == 0);
    
    fprintf(stderr, "Creating %d sites...\n", nsites);
    delaunay_site_t *st = makesites(nsites, normal);

    fprintf(stderr, "Building delaunay...\n");
    quad_arc_t e = delaunay_build (st, nsites);
    
    fprintf(stderr, "Plotting delaunay...\n");
    
    char *prefix = "out/delrandom";
    char *tag = distr;
    int32_t capLines = 2;
    epswr_figure_t *eps = delaunay_plot_new_figure(prefix, tag, capLines, st, nsites); 
    
    epswr_text(eps, "Delaunay and Voronoi", FALSE, 0.5, TRUE, FALSE);
    char *capt2 = NULL;
    asprintf(&capt2, "%s random sites", distr);
    epswr_text(eps, capt2, FALSE, 0.5, TRUE, FALSE);
    free(capt2);

    delaunay_plot_diagram(eps, e, st, nsites);
    epswr_end_figure(eps);
    return(0);
  }

delaunay_site_t *makesites(int32_t nsites, bool_t normal)
  { delaunay_site_t *st = notnull(malloc(nsites*sizeof(delaunay_site_t)), "no mem");
    
    srandom(4615);

    int32_t k;
    for (k=0; k < nsites; k++) 
      { st[k].index = k;
        /* Generate a random point {p} in the square {[-1 _ =1]^2}: */
        r2_t p;
        if (normal) 
          { /* Central limit approx to normal distribuition: */
            double s = 0, t = 0;
            int32_t j;
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
            deldebug_print_site ("given", &(st[k]));  fprintf(stderr, "\n");
            assert(abs(st[k].pt.c.c[0]) <= MC);
            assert(abs(st[k].pt.c.c[1]) <= MC);
            assert(abs(st[k].pt.c.c[2]) <= MC);
          }
      }
    return st;
  }
