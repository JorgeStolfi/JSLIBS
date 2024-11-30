#define PROG_NAME "test_r2_align_multiscale"
#define PROG_DESC "test of {r2_align_multiscale.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-08 11:16:27 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_align_COPYRIGHT \
  "Copyright � 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <assert.h>
#include <affirm.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <bool.h>
#include <r2x2.h>
#include <r2.h>
#include <ix.h>
#include <wt_table.h>
#include <wt_table_hann.h>

#include <r2_align.h>

#include <r2_align_multiscale.h>

int32_t main(int32_t argn, char **argv);

typedef double procedural_image_t(int32_t i, int32_t scale, double x, double y); 
  /* Evaluates image number {i} reduced to scale {scale} at point {x,y}. */

void test_align_one
  ( char *fname, 
    bool_t bias, 
    bool_t bal,
    bool_t monoscale, 
    bool_t quadopt, 
    int32_t NI, 
    r2_t cmp_rad, 
    r2_t adj_rad
  );
  /* Tests alignment algorithms for {NI} images derived by mapping of the
    mother function {fname}.  Uses window radius {cmp_rad},
    search radius {adj_rad}, balanced option {bal}. */

void test_align_choose_optimum(int32_t NI, r2_t pini[], r2_t adj_rad, r2_t popt[]);
  /* Stores into {popt[0..NI-1]} the optimum solution for the test
    with {NI} images, initial guess {pini[0..NI-1]}, and search radius {adj_rad}.
    makes sure that {popt-pini} is balanced. */

bool_t test_align_check_result
  ( int32_t NI, 
    r2_t psol[], 
    double Q2sol,
    r2_t pini[], 
    r2_t popt[], 
    r2_t rad[],
    double tol
  ); 
  /* Check correctness of solution {psol} against known optimum {popt} (if not NULL).
    The solution is correct if every point {psol[i]} differs from {popt[i]}
    by the same delta vector. Also checks that the delta vectors from {pini[i]}
    to {psol[i]} are within the stated search radii. */

void test_image_debug_points
  ( char *label, 
    int32_t NI, 
    r2_t p[], 
    int32_t scale, 
    r2_t pini[], 
    r2_t popt[], 
    r2_t rad[],
    double f2p, 
    double b2p
  );
  /* Prints the points {p[0..NI-1]} and the corresponding raw goal function value {f2p}
    and bias term {b2p}.  If {pini} is not NULL, prints also the difference between 
    each {p[i]} and the corresponding {pini[i]}.  Ditto if {popt} is not NULL */

void test_image_plot_goal
  ( r2_align_multiscale_mismatch_t *f2, 
    int32_t NI, 
    int32_t scale, 
    int32_t hwx, 
    double wx[], 
    int32_t hwy, 
    double wy[], 
    r2_t p[],
    r2_t *adj_rad
  );
  /* Writes a  file "out/f2.dat" with a random 2D slice of the goal function
    {f2(NI, scale, hwx, wx, hwy, wy, q)}, where {q} ranges in the neighborhood
    of {p}, with each {q[i]} inside the rectangle defined by {adj_rad}. */

double test_align_compute_image_mismatch_sqr
  ( int32_t ni,                 /* Number of images being compared. */
    procedural_image_t *eval,   /* procedural image evaluator. */
    int32_t scale,              /* Reduction scale. */
    int32_t hwx,                /* Half-width of comparison window. */
    double dx,                  /* X sampling step. */
    double wx[],                /* Horizontal weight table. */
    int32_t hwy,                /* Half-height of comparison window. */
    double dy,                  /* Y sampling step. */
    double wy[],                /* Vertical weight table. */
    r2_t p[]                    /* Sampling grid center for each image. */
  );
  /* Computes the mismatch between {ni} procedural images in the neighborhood of 
    points {p[0..ni-1]}.  
    
    Assumes that the value of image number {i} is given by
    [eval(i,x,y)}. Each image {i} is evaluated in a grid of sample
    points with {nx = 2*hwx+1} columns and {ny = 2*hwy+1} rows
    centered at the point {p[i]}. 
    
    For each position {ix,iy} in this grid, the variance of the
    sampled values is computed. The result is the average of those
    variances over the entire grid. */

void test_align_compute_avg_var(int32_t ni, double u[], double *avgP, double *varP);
  /* Computes the mean {*avgP} and variance {*varP} of the samples
    {u[0..ni-1]}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615*417);
    
    /* Parse arguments: */
    char *rest;

    int32_t na = 1;                                      /* Current argument index. */
    
    char *fname = argv[na++];                        /* Function name. */

    int32_t NI = (int32_t)strtol(argv[na++],&rest,10);            /* Number of images to align. */
    fprintf(stderr, "NI = %d\n", NI);
    demand((*rest) == 0, "invalid NI");

    bool_t bias = (int32_t)strtol(argv[na++],&rest,10);       /* 1 to add squared distance bias. */
    fprintf(stderr, "bias = %d\n", bias);
    demand((*rest) == 0, "invalid bias");

    bool_t monoscale = (int32_t)strtol(argv[na++],&rest,10);  /* 1 to perform single-scale enumeration. */
    fprintf(stderr, "monoscale = %d\n", monoscale);
    demand((*rest) == 0, "invalid monoscale");

    bool_t quadopt = (int32_t)strtol(argv[na++],&rest,10);    /* 1 to use quadratic optimization. */
    fprintf(stderr, "quadopt = %d\n", quadopt);
    demand((*rest) == 0, "invalid quadopt");

    r2_t cmp_rad;
    cmp_rad.c[0] = strtod(argv[na++],&rest);         /* Radius of comparison window. */
    cmp_rad.c[1] = strtod(argv[na++],&rest);         /* Radius of comparison window. */
    demand((*rest) == 0, "invalid window radius");

    r2_t adj_rad;
    adj_rad.c[0] = strtod(argv[na++],&rest);         /* Search radius. */
    adj_rad.c[1] = strtod(argv[na++],&rest);         /* Search radius. */
    demand((*rest) == 0, "invalid search radius");

    demand(argc == na, "wrong argument count");
    
    bool_t bal = TRUE;
    test_align_one(fname, bias, bal, monoscale, quadopt, NI, cmp_rad, adj_rad);
    return 0;
  }
  
void test_align_one
  ( char *fname, 
    bool_t bias, 
    bool_t bal,
    bool_t monoscale, 
    bool_t quadopt, 
    int32_t NI, 
    r2_t cmp_rad, 
    r2_t adj_rad
  )
  { 
    fprintf(stderr, "testing with %d procedural images...\n", NI);
    fprintf(stderr, "function = %s  bias = %c bal = %c\n", fname, "FT"[bias], "FT"[bal]);
    
    bool_t debug_points = FALSE; /* TRUE to print every probe point. */
    
    r2_t rad[NI];   /* Search radius for each coordinate. */
    double tol;     /* Search tolerance for each coordinate at each scale. */
    
    r2_t popt[NI];  /* Actual optimum alignment. */
    r2_t pini[NI];  /* Initial alignment. */

    /* Create the window weight tables: */
    int32_t hwx = (int32_t)ceil(fabs(cmp_rad.c[0])), nwx = 2*hwx + 1;
    double wx[nwx];
    wt_table_hann_fill(nwx, 0.0, wx, NULL);
    wt_table_normalize_sum(nwx, wx);
    
    int32_t hwy = (int32_t)ceil(fabs(cmp_rad.c[1])), nwy = 2*hwy + 1;
    double wy[nwy];
    wt_table_hann_fill(nwy, 0.0, wy, NULL);
    wt_table_normalize_sum(nwy, wy);
    
    int32_t NX = 20, NY = 10;   /* Nominal domain dimensions. */
    bool_t cmp_opt;         /* Should the solution be compared to {popt}? */
    
    double bias_factor = 1.0e-7; /* Scale factor for initial-point bias. */
    
    auto double f2_indiff(int32_t ni, r2_t p[], i2_t iscale);
      /* A goal function that returns 1.0 always.  If {bias} is true,
        the minimum of this function should be at {pini}. */
    
    auto double f2_images(int32_t ni, r2_t p[], i2_t iscale);
      /* A goal function that compares actual procedural images {img[0..NI-1][scale]},
        after aligning the points {p[0..NI-1]}. The images are compared
        within a weighted window of size {cmp_rad} centered on the aligned points.
        Requires that {iscale} is the same {scale} along both axes.
        If {bias} is FALSE, the minimum of this function should be near {popt}. */
    
    auto double f2_points(int32_t ni, r2_t p[], i2_t iscale);
      /* A goal function that ignores the images and simply returns
        the sum of the squared differences from {p[0..NI-1]} to the
        optimum {popt[0..NI-1]} divided by the respective radii.
        Coordinates with null {rad[i]} are skipped.
        
        If {bias} is FALSE, the minimum of this function
        should be very close to {popt}. */
    
    r2_align_multiscale_mismatch_t *f2_raw = NULL; /* The raw goal function, withot bias. */
    
    auto double Q2(int32_t ni, r2_t p[], i2_t iscale);
      /* The goal function that evaluates {f2_raw(ni, scale, p)} and
        adds the bias term if {bias} is true.  Also debugs the point 
        if {debug_points} is true. */
    
    auto double img_eval(int32_t i, int32_t scale, double x, double y);
      /* A goal function that evaluates one sample of procedural image number {i}
        reduced to the given {scale} at the point {x,y}. */
    
    /* The mother image: */
    /* Data for constituent waves: */
    int32_t NT = 10;      /* Number of waves. */
    double amp[NT];   /* Raw amplitude. */
    double phi[NT];   /* Initial wave phase (fraction of cycle). */
    r2_t frq[NT];     /* Wave frequency in X and Y (cycles per pixel). */

    auto double mother_image(int32_t scale, double x, double y);
      /* A procedural image which is the sum of several sine waves
        with various wavelengths, reduced to the given {scale}.
        Assumes the tables {NT,amp,phi,frq} have been initialized. */
    
    auto void initialize_mother_image(void);
      /* Initalizes {NT,amp,phi,frq} for {mother_image}. */
    
    /* Procedure body: */
    initialize_mother_image();
       
    /* Choose the goal function {f2_raw} according to {fname}. */
    /* Also set auxiliary parameters {pini[],popt[]}: */
    /* Choose the initial guess {pini}, radii {rad}, tolerance {tol}: */
    tol = 0.25;
    for (uint32_t i = 0;  i < NI; i++) 
      { pini[i] = (r2_t){{ 0.5*NX, 0.5*NY }};
        rad[i].c[0] = (i == 1 ? 0.0 : 0.5*(1 + sin(i))*adj_rad.c[0]);
        rad[i].c[1] = 0.5*(1 + cos(i))*adj_rad.c[1];
      }
        
    if (strcmp(fname, "indiff") == 0)
      { f2_raw = &f2_indiff;
        /* For this function, the optimum {popt} is the initial point {pini}: */
        for (uint32_t i = 0;  i < NI; i++) { popt[i] = pini[i]; }
        cmp_opt = bias; /* If {bias} is true, the solution should be {popt}. */
      }
    else if (strcmp(fname, "points") == 0)
      { f2_raw = &f2_points;
        /* Choose an optimum {popt} near the initial point {pini}: */
        test_align_choose_optimum(NI, pini, adj_rad, popt);
        cmp_opt = ! bias; /* If {bias} is false, the solution should be {popt}. */
     }
    else if (strcmp(fname, "images") == 0)
      { f2_raw = &f2_images;
        /* Choose the optimum alignment {popt} near the initial point {pini}: */
        test_align_choose_optimum(NI, pini, adj_rad, popt);
        cmp_opt = ! bias; /* If {bias} is false, the solution should be {popt}. */
      }
    else 
      { demand(FALSE, "invalid function name"); }
    
    /* Print raw function value and bias term for initial point: */
    i2_t iscale0 = (i2_t){{ 0, 0 }};
    double f2ini = f2_raw(NI, pini, iscale0);
    double b2ini = (bias ? bias_factor * r2_align_rel_disp_sqr(NI, pini, pini, rad) : 0);
    test_image_debug_points("initial guess ", NI, pini, 0, NULL, NULL, rad, f2ini, b2ini);
    
    /* Print raw function value and bias term for optimum point: */
    double f2opt = f2_raw(NI, popt, iscale0);
    double b2opt = (bias ? bias_factor * r2_align_rel_disp_sqr(NI, popt, pini, rad) : 0);
    test_image_debug_points("actual optimum", NI, popt, 0, pini, NULL, rad, f2opt, b2opt);
    
    /* Plot the goal function in the neighborhood of the initial guess: */ 
    fprintf(stderr, "plotting goal function...\n");
    test_image_plot_goal(Q2, NI, 0, hwx, wx, hwy, wy, pini, rad);
    
    /* Call the optimizer: */
    fprintf(stderr, "optimizing");
    debug_points = TRUE;
    r2_t psol[NI];  /* Computed alignment vector. */
    double Q2sol;   /* Goal function at {psol} including bias. */
    for (uint32_t i = 0;  i < NI; i++) { psol[i] = pini[i]; }
    if (monoscale)
      { if (quadopt)
          { fprintf(stderr, " (single_scale_quadopt)...\n");
            r2_align_multiscale_single_scale_quadopt(NI, iscale0, Q2, rad, bal, tol, psol, &Q2sol);
          }
        else
          { fprintf(stderr, " (single_scale_enum)...\n");
            r2_align_multiscale_single_scale_enum(NI, iscale0, Q2, rad, tol, psol, &Q2sol);
          }
      }
    else
      { fprintf(stderr, " (multi_scale)...\n");
        r2_align_multiscale(NI, Q2, quadopt, rad, tol, psol, &Q2sol);
      }
    fprintf(stderr, "done optimizing.\n");
    
    /* Print raw function value and bias term for computed optimum: */
    double f2sol = f2_raw(NI, psol, iscale0);
    double b2sol = (bias ? bias_factor * r2_align_rel_disp_sqr(NI, psol, pini, rad) : 0);
    test_image_debug_points("computed optimum", NI, psol, 0, pini, popt, rad, f2sol, b2sol);
    assert(Q2sol == f2sol + b2sol);
    
    bool_t ok = test_align_check_result(NI, psol, Q2sol, pini, (cmp_opt ? popt : NULL), rad, tol);

    if(! ok) { exit(1); }
    
    return;
    
    double Q2(int32_t ni, r2_t p[], i2_t iscale)
      { demand(ni == NI, "duh?");
        assert(iscale.c[0] == iscale.c[1]);
        int32_t scale = iscale.c[0];
        double f2p = f2_raw(ni, p, iscale);
        double b2p = 0.0;
        if (bias)
          { /* Add a tiny term that pulls the function towards the initial point: */
            double d2ini = r2_align_rel_disp_sqr(ni, p,  pini, rad); 
            b2p = bias_factor * d2ini;
          }
        if (debug_points) 
          { test_image_debug_points("probe point   ", ni, p, scale, pini, popt, rad, f2p, b2p); }
        return f2p + b2p;
      }
    
    double f2_indiff(int32_t ni, r2_t p[], i2_t iscale)
      { double f2p = 1.0;
        return f2p;
      }
    
    double f2_points(int32_t ni, r2_t p[], i2_t iscale)
      { /* Compute sum of relative square distances from actual optimum {opt}: */
        double f2p = r2_align_rel_disp_sqr(ni, p, popt, rad);
        return f2p;
      }
      
    double f2_images(int32_t ni, r2_t p[], i2_t iscale)
      { double dx = 1.0;  /* Sampling step in {X}. */
        double dy = 1.0;  /* Sampling step in {Y}. */
        assert(iscale.c[0] == iscale.c[1]);
        int32_t scale = iscale.c[0];
        double f2p = test_align_compute_image_mismatch_sqr
          ( NI, &img_eval, scale, hwx, dx, wx, hwy, dy, wy, p );
        return f2p;
      }
    
    double img_eval(int32_t i, int32_t scale, double x, double y)
      {
        /* Returns the mother image {img} displaced so that the point {pini[i]}
          is at {popt[i]} and reduced to the given scale. */

        double fs = pow(0.5, scale);
        double dx = fs*(popt[i].c[0] - pini[i].c[0]);
        double dy = fs*(popt[i].c[1] - pini[i].c[1]);
        return mother_image(scale, x - dx, y - dy);
      }
    
    double mother_image(int32_t scale, double x, double y)
      {
        int32_t t;
        double sf = pow(0.5, scale);  /* Image scaling factor. */
        double f0 = 1.0*sf;           /* Cutoff frequency (cycles per pixel). */
        double sum = 0.0;
        for (t = 0; t < 10; t++)
          { /* Choose the frequency vector of term {t}: */
            double fx = frq[t].c[0];                /* X frequency (cycles per pixel). */
            double fy = frq[t].c[1];                /* Y frequency (cycles per pixel). */
            double tph = (fx*x + fy*y)/sf + phi[t]; /* Wave phase at {x,y} (as fraction of cycle). */
            double val = amp[t]*sin(2*M_PI*tph); /* Raw wave value. */
            /* Attenuate for the sampling kernel at current scale: */
            double att = exp(-0.5*(fx*fx + fy*fy)/f0*f0);
            sum += att*val;
          }
        return sum;
      }
      
    void initialize_mother_image(void) 
      {
        int32_t t;
        for (t = 0; t < NT; t++)
          { /* Choose the frequency vector of term {t}: */
            double a = t;                     /* Azimuth of wave direction (radians). */ 
            double r0 = 2.0*pow(M_SQRT2,t);   /* Wavelength at scale 0 (pixels). */
            double fx = cos(a)/r0;            /* X frequency at scale 0 (cycles per pixel). */
            double fy = sin(a)/r0;            /* Y frequency at scale 0 (cycles per pixel). */
            frq[t] = (r2_t){{ fx, fy }};      /* Frequency vector at scale 0 (cycles per pixel). */
            phi[t] = sin(t);                  /* Wave initial phase (as fraction of cycle). */
            amp[t] = cos(3.5*t);              /* Raw amplitude. */
          }
      }
  }
    
void test_align_choose_optimum(int32_t NI, r2_t pini[], r2_t adj_rad, r2_t popt[])
  {
    if (NI == 0) { return; }
    
    int32_t a;
    for (a = 0; a < 2; a++)
      { double ra = 0.50*adj_rad.c[a]; /* Search radius along axis {a}. */
        double da[NI];                 /* Adjustments along axis {a}. */
        /* Set {da[0..NI-1]} to random displcements: */
        double sumd = 0; /* Sum of all {da[0..NI-1]}. */
        for (uint32_t i = 0;  i < NI; i++) { da[i] = (2*drandom()-1)*ra; sumd += da[i]; }
        /* Adjust {da} so that it is balanced: */
        double corr = -sumd/NI;
        for (uint32_t i = 0;  i < NI; i++) { da[i] += corr; }
        /* Convert delta vector {da[0..NI-1]} to alignment vector {popt[0..NI-1].c[a]}: */
        for (uint32_t i = 0;  i < NI; i++) { popt[i].c[a] = pini[i].c[a] + da[i]; }
      }
  }
    
void test_image_plot_goal
  ( r2_align_multiscale_mismatch_t *f2, 
    int32_t NI, 
    int32_t scale, 
    int32_t hwx, 
    double wx[], 
    int32_t hwy, 
    double wy[], 
    r2_t p[],
    r2_t *adj_rad
  )
  {
    /* Choose unit vectors {u[0..NI-1]} in {R^2}, with sum {(0,0)}: */
    r2_t u[NI];
    /* Throw a bunch of random unit vectors: */
    r2_t usum = (r2_t){{ 0, 0 }};
    for (uint32_t i = 0;  i < NI; i++) { r2_throw_dir(&(u[i])); r2_add(&(u[i]), &usum, &usum); }
    while (r2_norm(&usum) > 1.0e-5)
      { /* Shift {u[0..NI-1]} so that it has zero sum, and normalize: */
        r2_t uavg;
        r2_scale(1.0/NI, &usum, &uavg);
        usum = (r2_t){{ 0, 0 }};
        for (uint32_t i = 0;  i < NI; i++)
          { r2_sub(&(u[i]), &uavg, &(u[i]));
            double u2 = r2_dir(&(u[i]), &(u[i]));
            if (u2 < 1.0e-4) { r2_throw_dir(&(u[i])); }
            r2_add(&(u[i]), &usum, &usum);
          }
      }
      
    /* Choose unit vectors {v[i]} in {R^2}, with sum {(0,0)}, perp to {u[i]}: */
    r2_t v[NI];
    for (uint32_t i = 0;  i < NI; i++) { r2_cross(&(u[i]), &(v[i])); }
    
    /* Sweep the {p,u,v} plane and plot: */
    char *fname = "out/f2.dat";
    FILE *fpl = open_write(fname, TRUE);
    
    int32_t ns = 20; /* Number of steps in each direction, in each sense. */
    r2_t q[NI]; /* Probe points. */
    fprintf(stderr, "\n");
    for (int32_t iu = -ns; iu <= +ns; iu++)
      { double du = ((double)iu)/((double)ns);
        fprintf(stderr, ".");
        for (int32_t iv = -ns; iv <= +ns; iv++)
          { double dv = ((double)iv)/((double)ns);
            /* Compute the probe points {q[0..NI-1]}. */
            for (uint32_t i = 0;  i < NI; i++)
              { r2_mix(du,&(u[i]), dv, &(v[i]), &(q[i])); 
                r2_add(&(p[i]), &(q[i]), &(q[i])); 
              }
            /* Evaluate the function and plot: */
            i2_t iscale = (i2_t){{ scale, scale }};
            double f2p = f2(NI, q, iscale);
            fprintf(fpl, "%+9.6f %+9.6f  %12.6f\n", du, dv, f2p); 
          }
        /* Blank line between scanlines, for {gnuplot}: */
        fprintf(fpl, "\n"); 
      }

    fprintf(stderr, "\n");
    fclose(fpl);
  }

bool_t test_align_check_result
  ( int32_t NI, 
    r2_t psol[], 
    double Q2sol,
    r2_t pini[], 
    r2_t popt[], 
    r2_t rad[],
    double tol
  )
  {
    auto void compare_to_point(char *tag, r2_t pref[]);
      /* Prints the RMS absolute difference {psol[i] - pref[i]}, and 
        the RMS difference {psol[i] - pref[i]} relative to {rad[i]},
        over {i} in {0..NI-1}.  Ignores coordinate {j} of {psol[i],pref[i]} iff
        {rad[i].c[j]} is zero. */
    
    int32_t ind = 0; /* Indentation */
    fprintf(stderr, "%*schecking the solution\n", ind, "");
    
    bool_t valid = TRUE;
    if (popt != NULL) { compare_to_point("optimum", popt); }
    if (pini != NULL) { compare_to_point("initial", pini); }
    
    fprintf(stderr, "\n");
    return valid;
    
    void compare_to_point(char *tag, r2_t pref[])
      { int32_t nv = 0;        /* Number of nonzero coordinates in {rad[0..ni-1]}. */
    
        double sum_d2 = 0; /* Total abs squared disp between {psol} and {pref}. */
        double sum_e2 = 0; /* Total rel squared disp between {psol} and {pref}. */

        for (uint32_t i = 0;  i < NI; i++)
          { r2_t *q = &(psol[i]);
            r2_t *o = &(pref[i]);
            r2_t *r = &(rad[i]);
            r2_t d, e;
            r2_sub(q, o, &d);
            int32_t j;
            for (j = 0; j < 2; j++)
              { if (r->c[j] == 0) 
                  { d.c[j] = 0; e.c[j] = 0; }
                else
                  { e.c[j] = d.c[j]/r->c[j]; nv++; }
              }
            sum_d2 += r2_norm_sqr(&d);
            sum_e2 += r2_norm_sqr(&e);
          }
        double rms_d2 = (nv == 0 ? 0.0 : sqrt(sum_d2/nv));
        double rms_e2 = (nv == 0 ? 0.0 : sqrt(sum_e2/nv));
        /* Check whether delta vector is balanced: */
        fprintf(stderr, "%*s  Abs and rel RMS delta vector from %s = %10.6f %10.6f\n", ind, "", tag, rms_d2, rms_e2);
      }
  }

void test_image_debug_points
  ( char *label, 
    int32_t NI, 
    r2_t p[], 
    int32_t scale, 
    r2_t pini[], 
    r2_t popt[], 
    r2_t rad[],
    double f2p, 
    double b2p
  )
  { 
    auto void show_disp(char *tag, r2_t *q, r2_t *o, r2_t *r);
      /* Shows the difference {q} to {o}, abslolute and relative. */

    int32_t ind = 2*scale; /* Indentation */
    double fscale = pow(2.0, scale);
    fprintf(stderr, "%*s%s\n", ind, "", label);
    for (uint32_t i = 0;  i < NI; i++)
      { fprintf(stderr, "%*s  p[%02d] = ( %9.4f %9.4f )*2^%d", ind, "", i, p[i].c[0], p[i].c[1], scale);
        r2_t q; r2_scale(fscale, &(p[i]), &q);
        fprintf(stderr, " = ( %9.4f %9.4f )", q.c[0], q.c[1]);
        if (pini != NULL) { show_disp("ini", &q, &(pini[i]), &(rad[i])); }
        if (popt != NULL) { show_disp("opt", &q, &(popt[i]), &(rad[i])); }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "%*s  f2(p) = %22.10f  b2(p) = %22.10f  Q2(p) = %22.10f\n", ind, "", f2p, b2p, f2p+b2p);
    fprintf(stderr, "\n");
    return;

    void show_disp(char *tag, r2_t *q, r2_t *o, r2_t *r) 
      {
        r2_t d; r2_sub(q, o, &d);
        fprintf(stderr, " = %s + ( %+9.4f  %+9.4f )", tag, d.c[0], d.c[1]);
        double rdx = (r->c[0] == 0 ? 0.0 : d.c[0]/r->c[0]);
        double rdy = (r->c[1] == 0 ? 0.0 : d.c[1]/r->c[1]);
        fprintf(stderr, " = %s + ( %+9.4f  %+9.4f )*rad", tag, rdx, rdy);
      }
   }

double test_align_compute_image_mismatch_sqr
  ( int32_t ni,                 /* Number of images being compared. */
    procedural_image_t *eval,   /* Procedural image evaluator. */
    int32_t scale,              /* Reduction scale. */
    int32_t hwx,                /* Half-width of comparison window. */
    double dx,                  /* X sampling step. */
    double wx[],                /* Horizontal weight table. */
    int32_t hwy,                /* Half-height of comparison window. */
    double dy,                  /* Y sampling step. */
    double wy[],                /* Vertical weight table. */
    r2_t p[]                    /* Sampling grid center for each image. */
  )  
  { 
    double sum_wht_var = 0;
    double sum_wht = 0;
    int32_t jx, jy;
    for (jy = -hwy; jy <= +hwy; jy++)
      { double vy = jy*dy;     /* Vert offset of sample point. */
        for (jx = -hwx; jx <= +hwx; jx++)
          { double vx = jx*dx;    /* Horiz offset of sample point. */
            /* Sample the images at this grid sampling point: */
            double u[ni];
            for (uint32_t i = 0;  i < ni; i++)
              { /* Get samples of image {i}: */
                u[i] = eval(i, scale, p[i].c[0]+vx, p[i].c[1]+vy);
              }

            /* Compute mean and variance, and sampling point weight {wht}: */
            double avg, var;
            test_align_compute_avg_var(ni, u, &avg, &var);
            double wht = wx[jx+hwx]*wy[jy+hwy];
            sum_wht_var += wht*var;
            sum_wht += wht;
          }
      }
    return (sum_wht == 0 ? 0.0 : sum_wht_var / sum_wht);
  }

void test_align_compute_avg_var(int32_t ni, double u[], double *avgP, double *varP)
  {
    /* Compute the average {avg}: */
    double sum_u = 0; /* Sum of {u[ki]} */
    for (uint32_t i = 0;  i < ni; i++) { sum_u += u[i]; }
    double avg = sum_u/ni; 
    /* Compute the sample variance: */
    double sum_du2 = 0;
    for (uint32_t i = 0;  i < ni; i++) { double dui = u[i] - avg; sum_du2 += dui*dui; }
    double var = sum_du2/ni;
    (*avgP) = avg;
    (*varP) = var;
  }
    
            
    
