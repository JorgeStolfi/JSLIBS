/* See {minn_plot.h} */
/* Last edited on 2023-03-27 18:13:38 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <vec.h>
#include <rn.h>
#include <rmxn.h>
#include <rmxn_extra.h>
#include <jsfile.h>
#include <affirm.h>
#include <float_image.h>

#include <minn_plot.h>

#define Pr fprintf
#define Er stderr

#define minn_plot_RAD_FUDGE (1.0e-10)
  /* Perturbation to include edges of domain in the point grid. */
  
typedef void minn_plot_2D_proc_t(int32_t i0, int32_t i1, double Fy, double y[]);
  /* Type of a procedure that uses one data point of a 2D plot.
    The sample indices are {i0,i1}, and {y[0..n-1]}. Note that 
    not all combinations of {i0,i1} are generated. */
    
int32_t minn_plot_num_samples(double rad, double step);
  /* Number of samples that span the range {[-rad _ +rad]} with
    increment {step}. */

void minn_plot_print_vector(FILE *wr, char *name, int32_t n, double u[], double rad);
  /* Prints the plot direction {u[0..n-1]} and the associated radius {rad} to {wr},
    with label {name}. */

void minn_plot_2D_gen
  ( int32_t n, 
    double org[], 
    double u0[],
    double rad0,
    double u1[],
    double rad1, 
    bool_t box,
    double step,
    minn_goal_t *F,
    minn_plot_2D_proc_t *use
  );
  /* Like {minn_plot_2D_gnuplot}, but instead of writing each data point to {wr},
    calls {use(i0,i1,Fy,y)}. */
/* IMPLEMENTATIONS */

int32_t minn_plot_num_samples(double rad, double step)
  { int32_t NS = (int32_t)floor(rad/step + minn_plot_RAD_FUDGE);
    return NS;
  }

void minn_plot_1D_gnuplot
  ( FILE *wr, 
    int32_t n, 
    double org[], 
    double rad,
    double step, 
    minn_goal_t *F
  )
  {
    bool_t debug = TRUE;
    demand(n >= 1, "invalid {n}");
    
    /* Choose the plot directions: */
    int32_t nu = 2*n-1; /* Number of directions. */
    assert(nu > 0);
    double U[nu*n]; /* The rows are the directions. */
    rmxn_throw_directions(nu, n, U);
    if (debug)
      { Pr(Er, "plot directions:\n");
        for (int32_t i = 0; i < nu; i++)
          { double *Ui = &(U[i*n]);
            char *uname = NULL;
            asprintf(&uname, "u[%03d]", i);
            minn_plot_print_vector(Er, uname, n, Ui, rad);
            free(uname);
          }
      }
    
    /* Generate plot data for each direction: */
    double y[n];
    Pr(wr, "# n = %d\n", n);
    Pr(wr, "# fields: k e F(x) x[0] x[1] ...x[%d]\n", n-1);
    /* Compute number of samples on each side of 0 along each direction: */
    int32_t NS = minn_plot_num_samples(rad, step);
    for (int32_t i = 0; i < nu; i++)
      { /* Write plot data for direction {i}: */
        double *Ui = &(U[i*n]);
        for (int32_t k = -NS; k <= NS; k++)
          { double e = k * step;
            Pr(wr, "%+4d %14.6e ", k, e); 
            if (org == NULL)
              { rn_zero(n, y); }
            else
              { rn_copy(n, org, y); }
            rn_mix_in(n, e, Ui, y);
            double Fy = F(n, y);
            Pr(wr, "  %14.6e ", Fy); 
            for (int32_t j = 0; j < n; j++)
              { Pr(wr, " %14.6e", y[j]); }
            Pr(wr, "\n");
          }
        /* Separate plots by a blank line: */
        Pr(wr, "\n");
      }
    fflush(wr);
  }

void minn_plot_2D_gen
  ( int32_t n, 
    double org[], 
    double u0[],
    double rad0,
    double u1[],
    double rad1, 
    bool_t box,
    double step,
    minn_goal_t *F,
    minn_plot_2D_proc_t *use
  )
  {
    bool_t debug = TRUE;
    demand(n >= 2, "invalid {n}");
    
    if (debug)
      { Pr(Er, "plot directions:\n");
        minn_plot_print_vector(Er, "u0", n, u0, rad0);
        minn_plot_print_vector(Er, "u1", n, u1, rad1);
      }
    
    /* Compute number of samples on each side of 0 along each direction: */
    int32_t NS0 = minn_plot_num_samples(rad0, step);
    int32_t NS1 = minn_plot_num_samples(rad1, step);

    double y[n];  /* Working parameter vector. */

    for (int32_t i1 = -NS1; i1 <= NS1; i1++)
      { double e1 = i1*step;
        for (int32_t i0 = -NS0; i0 <= +NS0; i0++)
          { double e0 = i0*step;
            /* Check if sample point is to be plotted: */
            double ok = TRUE;
            if (! box)
              { /* Domain is an ellipse: */
                double d0 = e0/rad0;
                double d1 = e1/rad1;
                if (hypot(d0, d1) > 1.0 + 2*minn_plot_RAD_FUDGE) { ok = FALSE; }
              }
            if (ok)
              { /* Output point: */
                if (org == NULL)
                  { rn_zero(n, y); }
                else
                  { rn_copy(n, org, y); }
                rn_mix_in(n, e0, u0, y);
                rn_mix_in(n, e1, u1, y);
                double Fy = F(n, y);
                use(i0, i1, Fy, y);
              }
          }
      }
  }

void minn_plot_2D_gnuplot
  ( FILE *wr, 
    int32_t n, 
    double org[], 
    double u0[],
    double rad0,
    double u1[],
    double rad1, 
    bool_t box,
    double step,
    minn_goal_t *F
  )
  {
    int32_t i1_last = INT32_MAX; /* Last value of {i1} seen by {use}. */
    
    Pr(wr, "# n = %d\n", n);
    Pr(wr, "# fields: i0 i1 e0 e1 F(x) x[0] x[1] ...x[%d]\n", n-1);

    auto void use(int32_t i0, int32_t i1, double Fy, double y[]);

    minn_plot_2D_gen(n, org, u0, rad0, u1, rad1, box, step, F, &use);
    fflush(wr);
    
    return;
    
    /* INTERNAL IMPLS */
    
    void use(int32_t i0, int32_t i1, double Fy, double y[])
      { if ((i1 != i1_last) && (i1_last != INT32_MAX))
          { /* Separator line for {gnuplot}'s {splot}: */
            Pr(wr, "\n");
          }
        double e0 = i0*step;
        double e1 = i1*step;
        Pr(wr, "%+4d%+4d  %14.6e %14.6e ", i0, i1, e0, e1); 
        Pr(wr, "  %14.6e ", Fy); 
        for (int32_t j = 0; j < n; j++)
          { Pr(wr, " %14.6e", y[j]); }
        Pr(wr, "\n"); 
        i1_last = i1;
      }

  }

float_image_t *minn_plot_2D_float_image
  ( int32_t n, 
    double org[], 
    double u0[],
    double rad0,
    double u1[],
    double rad1, 
    bool_t box,
    double step,
    minn_goal_t *F
  )
  {
    /* Compute number of samples on each side of 0 along each direction: */
    int32_t NS0 = minn_plot_num_samples(rad0, step);
    int32_t NS1 = minn_plot_num_samples(rad1, step);

    int32_t NC = 1;
    int32_t NX = 2*NS0+1;
    int32_t NY = 2*NS1+1;
    float_image_t *img = float_image_new(NC, NX, NY);
    float_image_fill_channel(img, 0, NAN);
    
    auto void use(int32_t i0, int32_t i1, double Fy, double y[]);

    minn_plot_2D_gen(n, org, u0, rad0, u1, rad1, box, step, F, &use);

    return img;
    
    /* INTERNAL IMPLS */
    
    void use(int32_t i0, int32_t i1, double Fy, double y[])
      { int32_t ix = NS0 + i0; assert((ix >= 0) && (ix < NX));
        int32_t iy = NS1 + i1; assert((iy >= 0) && (iy < NY));
        float_image_set_sample(img, 0, ix, iy, (float)Fy);
      }
  }

void minn_plot_print_vector(FILE *wr, char *name, int32_t n, double u[], double rad)
  { fprintf(wr, "  %s = ", name);
    rn_print(wr, n, u);
    fprintf(wr, "  rad = %12.7f\n", rad);
  }
    
  
