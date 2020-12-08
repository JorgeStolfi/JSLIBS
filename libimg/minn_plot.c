/* See {minn_plot.h} */
/* Last edited on 2013-10-20 23:46:15 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <vec.h>
#include <rn.h>
#include <jsfile.h>
#include <float_image.h>

#include <minn_plot.h>

#define Pr fprintf
#define Er stderr

void minn_plot_2D_float_image
  ( FILE *wr, 
    int n,
    minn_plot_goal_t *F,
    double xo[],
    double xa[],
    double xb[],
    int NS
  )
  {
    int NC = 1;
    int NX = 2*NS+1;
    int NY = 2*NS+1;
    float_image_t *IMG = float_image_new(NC, NX, NY);
    
    /* allocate the arg vector {x}: */
    double *x = double_vec_new(n).e;
    
    /* Fill the plot image, remember the min point {EMin,HMin}: */
    int ia, ib;
    for (ia = -NS; ia <= +NS; ia++)
      { for (ib = -NS; ib <= +NS; ib++)
          { 
            double sa = ((double)ia)/NS; /* Relative abscissa of plot. */
            double sb = ((double)ib)/NS; /* Relative ordinate of plot. */
            rn_mix(n, sa, xa, sb, xb, x);
            rn_mix_in(n, 1-sa-sb, xo, x);
            double Fx = F(n, x);
            float_image_set_sample(IMG, 0, ia+NS, ib+NS, (float)Fx);
          }
      }
    
    float_image_write(wr, IMG);
    fflush(wr);
    float_image_free(IMG);
    free(x);
  }

void minn_plot_2D_gnuplot(FILE *wr, int n, minn_plot_goal_t *F, double x0[], int NS, double R)
  {
    /* Choose two orthogonal directions {u,v} in parameter space: */
    double u[n], v[n];
    if (n == 1)
      { u[0] = 1.0; v[0] = 0.0; }
    else if (n == 2)
      { u[0] = v[1] = 1.0; u[1] = v[0] = 0.0; }
    else if (n >= 3)
      { /* Pick a random direction for {u}: */
        rn_throw_dir(n, u);
        /* Pick a perpendicular direction for {v}: */
        double mv = 0;
        while (mv < 0.1)
          { rn_throw_dir(n, v);
            (void)rn_decomp(n, v, u, NULL, v);
            mv = rn_norm(n, v);
          }
        rn_dir(n, v, v);
      }
    else
      { assert(FALSE); }
        
    double xx[n];  /* Working parameter vector. */

    int Nu = (n < 1 ? 0 : NS);
    int Nv = (n < 2 ? 0 : NS);
    int iu, iv, k;
    fprintf(wr, "# n = %d\n", n);
    fprintf(wr, "# fields: fu fv F(x) x[0] x[1] ...x[%d]\n", n-1);
    for (iv = -Nv; iv <= Nv; iv++)
      { double fv = (Nv == 0 ? 0 : (iv*R)/Nv);
        for (iu = -Nu; iu <= +Nu; iu++)
          { double fu = (Nu == 0 ? 0 : (iu*R)/Nu);
            for (k = 0; k < n; k++) { xx[k] = x0[k] + fu*u[k] + fv*v[k]; }
            double Fx = F(n, xx);
            fprintf(wr, " %12.8f %12.8f ", fu, fv);
            fprintf(wr, " %+22.16e", Fx);
            for (k = 0; k < n; k++) { fprintf(wr, " %12.8f", xx[k]); }
            fprintf(wr, "\n");
          }
        if (Nv > 0) { fprintf(wr, "\n"); }
      }
    fflush(wr);
  }

void minn_plot_1D_gnuplot(FILE *wr, int n, minn_plot_goal_t *F, double x0[], int NS, double r)
  {
    double y[n];
    fprintf(wr, "# n = %d\n", n);
    fprintf(wr, "# fields: k e F(x) x[0] x[1] ...x[%d]\n", n-1);
    int k;
    for (k = -NS; k <= NS; k++)
      { double e = r*((double)k)/((double)NS);
        fprintf(wr, "%+4d %14.6e ", k, e); 
        rn_copy(n, x0, y);
        int i;
        for (i = 0; i < n; i++)
          { double sv = y[i];
            y[i] = x0[i] + e;
            double fy = F(n, y);
            fprintf(wr, " %14.6e", fy); 
            y[i] = sv;
          }
        fprintf(wr, "\n");
      }
  }
