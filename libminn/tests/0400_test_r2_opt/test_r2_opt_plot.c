/* See {test_r2_opt_plot.h}. */
/* Last edited on 2024-12-05 10:35:01 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsrandom.h>

#include <test_r2_opt_basic.h>

#include <test_r2_opt_plot.h>

void tr2o_choose_plot_directions(int32_t NI, r2_t rad[], r2_t u[], r2_t v[]);
  /* Chooses displacement vectors {u[0..NI-1],v[0..NI-1]} that define
    the grid of sampling points for {tr2o_plot_goal}.  Namely, sample point
    {q[ku,kv][i]} will be {ctr[i] + ku/NS*u[i] + kv/NS*v[i]}. */

double tr2o_throw_nonzero_U(void);
  /* Returns a number in {[-1 _ +1]} that is 
    definitely not zero. */

void tr2o_plot_goal
  ( r2_opt_goal_func_t *f2, 
    int32_t NI, 
    i2_t iscale, 
    r2_t ctr[],
    char *ctrtag,
    r2_t rad[]
  )
  {
    r2_t u[NI];
    r2_t v[NI];
    tr2o_choose_plot_directions(NI, rad, u, v);

    /* Sweep the {ctr,u,v} plane and plot: */
    char *fname = jsprintf("out/f2-%02d-%02d-%s.dat", iscale.c[0], iscale.c[1], ctrtag);
    FILE *fpl = open_write(fname, TRUE);
    
    int32_t NS = 20; /* Number of steps in each direction, in each sense. */
    r2_t q[NI]; /* Probe points. */
    fprintf(stderr, "\n");
    for (int32_t ku = -NS; ku <= +NS; ku++)
      { double du = ((double)ku)/((double)NS);
        fprintf(stderr, ".");
        for (int32_t kv = -NS; kv <= +NS; kv++)
          { double dv = ((double)kv)/((double)NS);
            /* Compute the probe points {q[0..NI-1]}. */
            for (uint32_t i = 0;  i < NI; i++)
              { r2_mix(du,&(u[i]), dv, &(v[i]), &(q[i])); 
                r2_add(&(ctr[i]), &(q[i]), &(q[i])); 
              }
            /* Evaluate the function and plot: */
            double f2p = f2(NI, q, iscale);
            fprintf(fpl, "%+16.10f %+16.10f  %22.16f\n", du, dv, f2p); 
          }
        /* Blank line between scanlines, for {gnuplot}: */
        fprintf(fpl, "\n"); 
      }

    fprintf(stderr, "\n");
    free(fname);
    fclose(fpl);
  }

void tr2o_choose_plot_directions(int32_t NI, r2_t rad[], r2_t u[], r2_t v[])
  { /* For each {i} in {0..NI-1}, we generate
      two orthogonal vectors {u[i],v[i]},
      that are nonzero along the variable coords of {ctr[i]}.
      If {p[i]} has no variable coords, both {u[i]} and {v[i]} 
      are zero.  If {p[i]} has only one variable coord, then 
      {u[i]} is nonzero along that coord, and {v[i]} is zero.
    */
    for (uint32_t i = 0;  i < NI; i++) 
      { /* Count variable coords {nv} of {p[i]}, find last nonzero var {jf}: */
        int32_t nv = 0;
        int32_t jf = -1;
        for (uint32_t j = 0;  j < 2; j++)
          { double rij = rad[i].c[j];
            if (rij > 0.0) { nv++; jf = j; }
          }
        /* Set {u[i]} and {v[i]} accordingly: */
        if (nv == 0)
          { /* Point {p[i]} is fixed, so {u[i]} and {v[i]} are both zero: */
            u[i] = (r2_t){{ 0.0, 0.0 }};
            v[i] = (r2_t){{ 0.0, 0.0 }};
          }
        else if (nv == 1)
          { /* Point {p[i]} can vary only along coord {jf}. */
            /* Set {u[i]} to a random multiple of {rad[i]}: */
            double rijf = rad[i].c[jf];
            u[i].c[jf] = rijf*tr2o_throw_nonzero_U();
            u[i].c[1-jf] = 0.0;
            v[i] = (r2_t){{ 0.0, 0.0 }};
          }
        else if (nv == 2)
          { /* Generate a random vector {ui} in {[-1 _ +1]^2}, nonzero along variable coords: */
            for (uint32_t j = 0;  j < 2; j++)
              { double rij = rad[i].c[j];
                u[i].c[j] = rij*tr2o_throw_nonzero_U();
              }
            /* Get the direction {dui} of that vector: */
            r2_t dui; 
            double mu = r2_dir(&(u[i]), &dui);
            assert(mu > 0.0);
            /* Get direction {dvi} orthogonal to that: */
            r2_t dvi = (r2_t){{ -dui.c[1], +dui.c[0] }};
            /* Compute the rough extent {m} of box {rad[i]} in direction {dvi}: */
            double me = rad[i].c[0]*fabs(dvi.c[0]) + rad[i].c[1]*fabs(dvi.c[1]);
            assert(me > 0.0);
            /* Set {v[i]} orthogonal to {u[i]} scaled by a random fraction of {me}: */
            double fi = tr2o_throw_nonzero_U();
            r2_scale(fi*me, &dvi, &(v[i]));
          }
        else
          { assert(FALSE); }
      }

    /* Make sure that at least one {u[i]} and one {v[i]} span the {rad[i]} box: */
    tr2o_plot_grid_normalize_to_span(NI, u, rad);
    tr2o_debug_points(0, "plot direction u", NI, "pltu", u, NULL, NULL, rad, NULL, NAN, NAN);
    tr2o_plot_grid_normalize_to_span(NI, v, rad);
    tr2o_debug_points(0, "plot direction v", NI, "pltv", v, NULL, NULL, rad, NULL, NAN, NAN);
  }

void tr2o_plot_grid_normalize_to_span(int32_t NI, r2_t u[], r2_t rad[])
  {
    double tmax = tr2o_compute_rel_span(NI, u, rad); 
    double sf = 1.0/tmax; /* Scaling factor. */
    for (uint32_t i = 0;  i < NI; i++) { r2_scale(sf, &(u[i]), &(u[i])); }
  }
  
double tr2o_compute_rel_span(int32_t NI, r2_t u[], r2_t rad[]) 
  { double tmax = 0.0;
    for (uint32_t i = 0;  i < NI; i++) 
      { /* Get the direction {dui} of {u[i]}: */
        r2_t dui; 
        double mu = r2_dir(&(u[i]), &dui);
        if (mu > 0.0)
          { /* Get elipse extent {me} along {dui}: */
            double me = rad[i].c[0]*fabs(dui.c[0]) + rad[i].c[1]*fabs(dui.c[1]);
            assert(me > 0.0);
            double ti = mu/me;
            if (ti > tmax) { tmax = ti; }
          }
      }
    return tmax;
  }

void tr2o_write_test_image
  ( tr2o_image_eval_proc_t *eval, 
    int32_t i,        /* Image index. */
    i2_t iscale,  /* Image shrink scale in each axis. */
    i2_t wsize,   /* Comparison widow size along each axis. */
    r2_t ctr,     /* Center of image (unscaled). */
    char *ctrtag, /* Type of center, e.g. "ini" or "opt". */
    r2_t rad      /* Search radius for this scale. */
  )
  {
    /* Get half-widths of sampling window: */
    demand((wsize.c[0] % 2) == 1, "window width must be odd");
    demand((wsize.c[1] % 2) == 1, "window height must be odd");
    int32_t hwx = (wsize.c[0]-1)/2; 
    int32_t hwy = (wsize.c[1]-1)/2;

    r2_t scale = (r2_t){{ pow(2.0, iscale.c[0]), pow(2.0, iscale.c[1]) }};
    
    char *fname = jsprintf("out/image-%02d-%02d-%03d-%s.pgm", iscale.c[0], iscale.c[1], i, ctrtag);
    /* Image size: */
    int32_t HX = 2*((int32_t)ceil(rad.c[0])+hwx)+1; int32_t NX = 2*HX + 1;
    int32_t HY = 2*((int32_t)ceil(rad.c[1])+hwy)+1; int32_t NY = 2*HY + 1;
    fprintf(stderr, "writing image file %s (%d x %d)\n", fname, NX, NY);
    
    /* We cannot import {libimg} yet so let's write ASCII PGM: */
    FILE *wr = open_write(fname, TRUE);
    int32_t maxval = 65535;
    fprintf(wr, "P2\n");
    fprintf(wr, "%d %d\n", NX, NY);
    fprintf(wr, "%d", maxval);
    for (int32_t iy = -HY; iy <= +HY; iy++)
      { double ysc = ctr.c[1]/scale.c[1] + iy;
        for (int32_t ix = -HX; ix <= +HX; ix++)
          { double xsc = ctr.c[0]/scale.c[0] + ix;
            double val = eval(i, iscale, xsc, ysc);
            assert((val >= 0.0) && (val <= 1.0));
            int32_t ival = (int32_t)floor(val*(maxval + 0.9999999));
            assert((ival >= 0) && (ival <= maxval));
            fprintf(wr, (((ix + HX) % 20) == 0 ? "\n" : " "));
            fprintf(wr, "%d", ival);
          }
          fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
    free(fname);
    fclose(wr);
  }

double tr2o_throw_nonzero_U(void)
  { double z;
    do { z = 2*drandom() - 1; } while (fabs(z) < 1.0e-4);
    return z;
  }

