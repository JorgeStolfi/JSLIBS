#define PROG_NAME "test_r2_opt"
#define PROG_DESC "test of {r2_opt.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-03-19 14:41:37 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define tr2o_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <r2_opt.h>
#include <wt_table.h>
#include <wt_table_hann.h>
#include <ix.h>
#include <r2.h>
#include <r2x2.h>
#include <bool.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <jsprintf.h>
#include <affirm.h>

#include <test_r2_opt_info.h>
#include <test_r2_opt_basic.h>
#include <test_r2_opt_plot.h>
#include <test_r2_opt_tools.h>

/* #include <float_image.h> */
/* #include <float_image_interpolate.h> */
/* #include <float_image_io_pnm.h> */

#define tr2o_MAX_SCALE (15)
  /* Hopefully the scale will never be higher than this. */ 

void tr2o_do_one_test
  ( char *fname, 
    bool_t bias, 
    bool_t monoscale, 
    bool_t quadopt, 
    uint32_t NI, 
    i2_t cmp_isc,
    r2_t adj_rad,
    r2_t adj_stp
  );
  /* Tests alignment algorithms for {NI} images derived by displacement of the
    mother image {fname}.  Uses 
    search radius {adj_rad}, and step/precision {adj_stp}. 
    
    If {monoscale} is {TRUE}, the the images are scaled by {1/2^cmp_isc.c[j]}
    in each axis {j}. If {monoscale} is {FALSE}, the {cmp_isc} 
    parameter should be {(0,0)}. */

void tr2o_make_weight_table(int32_t hw, int32_t *nwp, double **wtp);
  /* Makes a weight table for sampling a univariate function over a 
    symmetric window with {2*hw+1} samples, with unit steps. */
    
void tr2o_choose_pini_arad_astp
  ( uint32_t NI, 
    r2_t adj_rad, 
    r2_t adj_stp, 
    r2_t pini[], 
    r2_t arad[], 
    r2_t astp[]
  );
  /* Chooses the initial guess {pini[i]}, the search radius {arad[i]},
    and the search step {astp[i]}, for {i=0..NI-1}.  
    
    The intial guess is at the center of the nominal image domain.  
    The search radius and search step are {adj_rad} and {adj_stp},
    perturbed by some random amount for each image; or zero for
    some coordinates and images. */
    
void tr2o_choose_optimum(uint32_t NI, r2_t pini[], r2_t arad[], r2_t popt[]);
  /* Stores into {popt[0..NI-1]} the optimum solution for the test
    with {NI} images near the initial guess {pini[0..NI-1]}, 
    and search radii {arad[0..NI-1]}. */

bool_t tr2o_ckeck_result
  ( uint32_t NI, 
    r2_t psol[], 
    double Q2sol,
    r2_t pini[], 
    r2_t popt[], 
    r2_t arad[],
    r2_t astp[]
  ); 
  /* Check correctness of solution {psol} against known optimum {popt}
    (if not NULL). The solution is correct if every point {psol[i]}
    differs from {popt[i]} by the same displacement vector, with error
    smaller than {astp[i]}. Also checks that the displacements from
    {pini[i]} to {psol[i]} are within the stated search radius
    {arad[i]}. */

double tr2o_compute_image_mismatch_sqr
  ( uint32_t NI,                       /* Number of images being compared. */
    tr2o_image_eval_proc_t *eval, /* Image evaluator. */
    r2_t p[],     /* Unscaled sampling grid center for each image. */
    i2_t iscale,  /* Image shrink scale in each axis. */
    i2_t wsize,   /* Comparison widow size along each axis. */
    double wtx[], /* Weight table for X displacement. */
    double wty[]  /* Weight table for Y displacement. */
  );
  /* Computes the mismatch between {NI} images in the neighborhood of points
    {p[0..NI-1]}.  The images are implicitly shrunk (with aliasing)
    by factor {1/2^iscale.c[j]} along each axis {j}. Ditto for the
    points {p[i]}. 
    
    Each scaled image {i} is evaluated in a grid of sample points
    {(xsc,ysc)}, with unit steps, with {wsize.c[0]} columns and {wsize.c[1]}
    rows centered at the point {p[i]} implicitly scaled by the same
    factors. Assumes that the value of the shrunk image number {i} at
    the scaled point {(xsc,ysc)} is given by {eval(i,iscale,xsc,ysc)}.
    
    For each position {ix,iy} in this grid, the variance of the
    sampled values is computed. The result is the average of those
    variances over the entire grid. */

void tr2o_compute_avg_var(uint32_t nz, double z[], double *avgP, double *varP);
  /* Computes the mean {*avgP} and variance {*varP} of the samples
    {z[0..nz-1]}. */

int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615*417);
    
    /* Parse arguments: */
    fprintf(stderr, "=== command-line arguments =========================\n");
    char *rest;

    uint32_t na = 1;                                      /* Current argument index. */
    
    char *fname = argv[na++];                        /* Name of goal function for optimization. */
    fprintf(stderr, "fname = %s\n", fname);
    
    bool_t bias = (bool_t)strtol(argv[na++],&rest,10);       /* 1 to add squared distance bias. */
    demand(((*rest) == 0) && (bias >= 0) && (bias <= 1), "invalid bias");
    fprintf(stderr, "bias = %d\n", bias);

    bool_t monoscale = (bool_t)strtol(argv[na++],&rest,10);  /* 1 to perform single-scale enumeration. */
    demand(((*rest) == 0) && (monoscale >= 0) && (monoscale <= 1), "invalid monoscale");
    fprintf(stderr, "monoscale = %d\n", monoscale);

    bool_t quadopt = (bool_t)strtol(argv[na++],&rest,10);    /* 1 to use quadratic optimization. */
    demand(((*rest) == 0) && (quadopt >= 0) && (quadopt <= 1), "invalid quadopt");
    fprintf(stderr, "quadopt = %d\n", quadopt);

    uint32_t NI = (uint32_t)strtol(argv[na++],&rest,10);       /* Number of images to align. */
    demand(((*rest) == 0) && (NI > 0), "invalid NI");
    fprintf(stderr, "number of images = %d\n", NI);

    r2_t adj_rad;
    adj_rad.c[0] = strtod(argv[na++],&rest);         /* X search radius. */
    demand(((*rest) == 0) && (adj_rad.c[0] > 0.0), "invalid X search radius");
    adj_rad.c[1] = strtod(argv[na++],&rest);         /* Y search radius. */
    demand(((*rest) == 0) && (adj_rad.c[1] > 0.0), "invalid Y search radius");
    fprintf(stderr, "search radius = %4.1f %4.1f\n", adj_rad.c[0], adj_rad.c[1]);

    r2_t adj_stp;
    adj_stp.c[0] = strtod(argv[na++],&rest);         /* X enum step or opt precision. */
    demand(((*rest) == 0) && (adj_stp.c[0] >= 0.0), "invalid X search step");
    adj_stp.c[1] = strtod(argv[na++],&rest);         /* Y enum step or opt precision. */
    demand(((*rest) == 0) && (adj_stp.c[1] >= 0.0), "invalid Y search step");
    demand((adj_stp.c[0] == 0.0) == (adj_stp.c[1] == 0.0), "invalid X,Y steps");
    fprintf(stderr, "enum search step = %4.1f %4.1f\n", adj_stp.c[0], adj_stp.c[1]);

    i2_t cmp_isc;
    cmp_isc.c[0] = (int32_t)strtol(argv[na++],&rest,10);         /* X domain scale. */
    demand(((*rest) == 0) && (cmp_isc.c[0] >= 0), "invalid X domain scale");
    cmp_isc.c[1] = (int32_t)strtol(argv[na++],&rest,10);         /* Y domain scale. */
    demand(((*rest) == 0) && (cmp_isc.c[1] >= 0), "invalid Y domain scale");
    fprintf(stderr, "comparison scale = %d %d\n", cmp_isc.c[0], cmp_isc.c[1]);

    demand(argc == na, "wrong argument count");
    fprintf(stderr, "---------------------------------------\n");
    
    tr2o_do_one_test(fname, bias, monoscale, quadopt, NI, cmp_isc, adj_rad, adj_stp);
    return 0;
  }
  
void tr2o_do_one_test
  ( char *fname, 
    bool_t bias, 
    bool_t monoscale, 
    bool_t quadopt, 
    uint32_t NI, 
    i2_t cmp_isc,
    r2_t adj_rad,
    r2_t adj_stp
  )
  { 
    fprintf(stderr, "testing with %d images...\n", NI);
    fprintf(stderr, "mother image = %s  bias = %c\n", fname, "FT"[bias]);
    
    bool_t debug_points = FALSE; /* TRUE to print every probe point. */
    bool_t debug_opt = TRUE; /* Passed to {r2_opt.h} routines. */
    
    bool_t compare_to_popt;      /* Should the solution be compared to {popt}? */
    
    double bias_mag = 1.0e-7; /* Magnitude of initial-point bias. */
    
    /* Select the scale for debugging function values etc: */
    i2_t iscale = ( monoscale ? cmp_isc : (i2_t){{ 0, 0 }} );
    
    /* The window sampling weight tables for each coordinate: */
    i2_t cmp_rad = (i2_t){{ 3, 3 }};
    i2_t wsize; 
    double *wtx; tr2o_make_weight_table(cmp_rad.c[0], &wsize.c[0], &wtx);
    double *wty; tr2o_make_weight_table(cmp_rad.c[0], &wsize.c[1], &wty);

    /* Goal functions for {r2_opt}: */

    r2_opt_goal_func_t *f2_raw = NULL; /* */
      /* The raw goal function for optimization, without bias.
        Selected by {fname}. */
    
    uint32_t nf2; /* Counts calls to {nf2}. */
    
    auto double f2_full(uint32_t ni, r2_t p[], i2_t iscale);
      /* The goal function that evaluates {f2_raw(ni,p,iscale)} and
        adds the bias term {bias_term(ni,p)} if {bias} is true.  
        Also debugs the point if {debug_points} is true. */
    
    auto double bias_term(r2_t p[]);
      /* The bias term (small multiple of squared distance from {pini}. */
    
    /* Alternatives for {f2_raw}: */

    auto double f2_indiff(uint32_t ni, r2_t p[], i2_t iscale);
      /* A goal function that returns 1.0 always.  If {bias} is true,
        the minimum of this function should be at {pini}. */
    
    auto double f2_imgmis(uint32_t ni, r2_t p[], i2_t iscale);
      /* A goal function that compares the images {0..NI-1},
        in the neighborhood of the points {p[0..NI-1]}.
        The domains of the images are implicitly scaled by
        {1/2^iscale.c[j]} (with antialiasing) along each axis {j}.  
        Ditto for the coordinates of the the points {p[i]}.
        
        Specifically, evaluates each image {i}, implicitly scaled, at
        points in a grid of {wsize.c[0]} by {wsize.c[1]} sampling points
        {q[kx,xy][i]} with unit steps centered at {p[i]}, also scaled.
        That grid covers a rectangle with size {arad.c[j]} along each
        axis {j}. Let {f[kx,xy][i]} be the values of the {ni} images at
        those points, and {var[kx,ky]} the variance of the values for
        the same grid indices {kx,ky}. The value of {f2_imgmis} is the
        average of those variances, weighted by the window weights.
        
        If {bias} is FALSE, the minimum of this function
        should be near {popt}. */
    
    auto double f2_optdst(uint32_t ni, r2_t p[], i2_t iscale);
      /* A goal function that ignores the images and simply returns
        the sum of the squared differences from {p[0..NI-1]} to the
        optimum {popt[0..NI-1]} divided by the respective radii.
        Coordinates with null {arad[i]} are skipped.
        If {bias} is FALSE, the minimum of this function
        should be very close to {popt}. */
    
    /* images for {f2_imgmis}: */
    
    auto double image_eval(uint32_t i, i2_t iscale, double xsc, double ysc);
      /* A goal function that evaluates one sample of image number {i}
        with its domain implicitly reduced by {1/2^iscale.c[j]} along
        each axis {j}, evaluated at point {(xsc,ysc)}, assumed to be
        already scaled. */
    
    /* The mother image: */
    
    /* Data for constituent waves: */
    uint32_t mom_NF = 30;      /* Number of waves. */
    double mom_amp[mom_NF];   /* Raw amplitude. */
    r2_t mom_phi[mom_NF];   /* Initial wave phase (fraction of cycle). */
    r2_t mom_frq[mom_NF];     /* Wave frequency in X and Y (cycles per pixel). */
    
    /* Set auxiliary parameters {pini[],popt[],arad[],astp[]}: */
    r2_t popt[NI];  /* Actual optimum alignment. */
    r2_t pini[NI];  /* Initial alignment. */
    r2_t arad[NI];  /* Search radius for each coordinate. */
    r2_t astp[NI];  /* Search step or precision for each coordinate. */ 
    tr2o_choose_pini_arad_astp(NI, adj_rad, adj_stp, pini, arad, astp);
    
    /* Choose the goal function {f2_raw} for optimization according to {fname}. */

    if (strcmp(fname, "indiff") == 0)
      { f2_raw = &f2_indiff;
        /* For this function, the optimum {popt} is the initial point {pini}: */
        for (uint32_t i = 0;  i < NI; i++) { popt[i] = pini[i]; }
        compare_to_popt = bias; /* If {bias} is true, the solution should be {popt}. */
      }
    else if (strcmp(fname, "optdst") == 0)
      { f2_raw = &f2_optdst;
        /* Choose an optimum {popt} near the initial point {pini}: */
        tr2o_choose_optimum(NI, pini, arad, popt);
        compare_to_popt = ! bias; /* If {bias} is false, the solution should be {popt}. */
     }
    else if (strcmp(fname, "imgmis") == 0)
      { f2_raw = &f2_imgmis;
        /* Initialize the parameters of the mother image: */
        tr2o_tools_initialize_mother_image(mom_NF, mom_frq, mom_phi, mom_amp, TRUE);
        /* Choose the optimum displacements {popt} near the initial point {pini}: */
        tr2o_choose_optimum(NI, pini, arad, popt);
        compare_to_popt = ! bias; /* If {bias} is false, the solution should be {popt}. */
        
        /* Write the images out for debuging: */
        for (int32_t isc = -1; isc < 3; isc++)
          { i2_t ipscale = ( isc < 0 ? iscale : (i2_t){{ isc, isc }} ); /* Plotting scale. */
            for (uint32_t i = 0;  i < NI; i++)
              { tr2o_write_test_image(image_eval, i, ipscale, cmp_rad, pini[i], "ini", arad[i]);
                tr2o_write_test_image(image_eval, i, ipscale, cmp_rad, popt[i], "opt", arad[i]);
              }
          }
      }
    else 
      { demand(FALSE, "invalid goal function name"); }
    
    /* Print raw function value and bias term for initial point: */
    double f2ini = f2_raw(NI, pini, iscale);
    double b2ini = (bias ? bias_term(pini) : 0);
    tr2o_debug_points(0, "initial guess ", NI, "pini", pini, NULL, NULL, arad, astp, f2ini, b2ini);
    
    /* Print raw function value and bias term for optimum point: */
    double f2opt = f2_raw(NI, popt, iscale);
    double b2opt = (bias ? bias_term(pini) : 0);
    tr2o_debug_points(0, "actual optimum", NI, "popt", popt, pini, NULL, arad, astp, f2opt, b2opt);
    
    /* Plot the goal function in the neighborhood of the initial guess and optimum point: */ 
    for (int32_t isc = -1; isc < 3; isc++)
      { i2_t ipscale = ( isc < 0 ? iscale : (i2_t){{ isc, isc }} ); /* Plotting scale. */
        fprintf(stderr, "plotting goal function at scale (%d,%d)...\n", ipscale.c[0], ipscale.c[1]);
        nf2 = 0;
        tr2o_plot_goal(f2_full, NI, ipscale, pini, "ini", arad);
        fprintf(stderr, "%d calls to {f2_full} for plotting.\n", nf2);
        nf2 = 0;
        tr2o_plot_goal(f2_full, NI, ipscale, popt, "opt", arad);
        fprintf(stderr, "%d calls to {f2_full} for plotting.\n", nf2);
      }
      
    /* If search step is zero, just plot, skip optimization: */
    if ((adj_stp.c[0] == 0.0) && (adj_stp.c[1] == 0.0))
      { fprintf(stderr, "no optimization: search step is zero\n");
        return;
      }
    
    /* Call the optimizer: */
    fprintf(stderr, "optimizing");
    
    debug_points = TRUE;
    r2_t psol[NI];  /* Computed alignment vector. */
    double f2sol;   /* Goal function at {psol} including bias. */
    nf2 = 0;
    for (uint32_t i = 0;  i < NI; i++) { psol[i] = pini[i]; }
    if (monoscale)
      { if (quadopt)
          { fprintf(stderr, " (single_scale_quadopt)...\n");
            r2_opt_single_scale_quadopt(NI, iscale, f2_full, arad, astp, psol, &f2sol, debug_opt);
          }
        else
          { fprintf(stderr, " (single_scale_enum)...\n");
            r2_opt_single_scale_enum(NI, iscale, f2_full, arad, astp, psol, &f2sol, debug_opt);
          }
      }
    else
      { fprintf(stderr, " (multi_scale)...\n");
        r2_opt_multi_scale(NI, f2_full, quadopt, arad, astp, psol, &f2sol, debug_opt);
      }
    fprintf(stderr, "done optimizing.\n");
    fprintf(stderr, "%d calls to {f2_full}.\n", nf2);
    
    /* Print raw function value and bias term for computed optimum: */
    double f2_raw_sol = f2_raw(NI, psol, iscale);
    double bias_sol = (bias ? bias_term(psol) : 0);
    tr2o_debug_points(0, "computed optimum", NI, "psol", psol, pini, popt, arad, astp, f2_raw_sol, bias_sol);
    assert(f2sol == f2_raw_sol + bias_sol);
    
    /* Print again raw function value and bias term for optimum point: */
    tr2o_debug_points(0, "actual optimum", NI, "popt", popt, pini, NULL, arad, astp, f2opt, b2opt);
    
    bool_t ok = tr2o_ckeck_result(NI, psol, f2sol, pini, (compare_to_popt ? popt : NULL), arad, astp);

    if(! ok) { exit(1); }
    
    return;
    
    /* INTERNAL IMPLEMENTATIONS */
    
    double bias_term(r2_t p[])
      { double d2ini = r2_opt_rel_disp_sqr(NI, p, pini, arad, astp); 
        double b2p = bias_mag * d2ini;
        return b2p;
      }
    
    double f2_full(uint32_t ni, r2_t p[], i2_t iscale)
      { demand(ni == NI, "duh?");
        nf2++; 
        double f2p = f2_raw(ni,p,iscale);
        double b2p = 0.0;
        if (bias) { b2p = bias_term(p); }
        if (debug_points) 
          { uint32_t indent = 2*(uint32_t)imax(iscale.c[0], iscale.c[1]); /* Indentation. */
            tr2o_debug_points(indent, "probe point   ", ni, "pcur", p, pini, popt, arad, astp, f2p, b2p);
          }
        return f2p + b2p;
      }
    
    double f2_indiff(uint32_t ni, r2_t p[], i2_t iscale)
      { double f2p = 1.0;
        return f2p;
      }
    
    double f2_optdst(uint32_t ni, r2_t p[], i2_t iscale)
      { /* Compute sum of relative square distances from actual optimum {opt}: */
        double f2p = r2_opt_rel_disp_sqr(ni, p, popt, arad, astp);
        return f2p;
      }
      
    double f2_imgmis(uint32_t ni, r2_t p[], i2_t iscale)
      { double f2p = tr2o_compute_image_mismatch_sqr
          ( NI, &image_eval, p, iscale, wsize, wtx, wty );
        return f2p;
      }
    
    double image_eval(uint32_t i, i2_t iscale, double xsc, double ysc)
      {
        /* Returns the mother image {img}, displaced
          by {popt[j]}, shrunk by {iscale.c[j]} along axis {j},
          then evaluated at {(xsc,ysc)}. */
        r2_t scale = (r2_t){{ pow(2.0, iscale.c[0]), pow(2.0, iscale.c[1]) }};
        double dxsc = popt[i].c[0]/scale.c[0];
        double dysc = popt[i].c[1]/scale.c[1];
        r2_t psc = (r2_t){{ xsc - dxsc, ysc - dysc }};
        return tr2o_tools_eval_mother_image(&psc, &scale, mom_NF, mom_frq, mom_phi, mom_amp);
      }
    
  }

void tr2o_make_weight_table(int32_t hw, int32_t *nwp, double **wtp)
  { demand (hw >= 0, "invalid sampling radius");
    uint32_t nw = (uint32_t)(2*hw + 1);
    double *wt = notnull(malloc(nw*sizeof(double)), "no mem");
    wt_table_hann_fill(nw, 0.0, wt, NULL);
    wt_table_normalize_sum(nw, wt);
    (*nwp) = (int32_t)nw;
    (*wtp) = wt;
  }
    
void tr2o_choose_pini_arad_astp
  ( uint32_t NI, 
    r2_t adj_rad, 
    r2_t adj_stp, 
    r2_t pini[], 
    r2_t arad[], 
    r2_t astp[]
  )   
  {
    /* Coordinates {kfix + r*mfix} for integer {r} will be fixed: */
    uint32_t kfix = 1;
    uint32_t mfix = 7;
    /* Define the parameters: */
    uint32_t k = 0; /* Counts coordinates. */
    for (uint32_t i = 0;  i < NI; i++) 
      { /* Choose the initial guess {pini}: */
        pini[i] = (r2_t){{ 50.0, 50.0 }};
        
        for (uint32_t j = 0;  j < 2; j++)
          { if (((k - kfix) % mfix == 0) || (adj_rad.c[j] == 0))
              { /* Do not optimize this coordinate: */
                arad[i].c[j] = 0.0;
                astp[i].c[j] = 0.0;
              }
            else
              { /* Choose the search radius {arad[i].c[j]} by wriggling {adj_rad.c[j]}: */
                arad[i].c[j] = (0.85 + 0.15*sin(k))*adj_rad.c[j];
                if (adj_stp.c[j] == 0.0)
                  { /* Plot but do not optimize this coordinate: */
                    astp[i].c[j] = 0.0;
                  }
                else
                  { /* Choose the enum step {astp[i].c[j]} by wriggling {adj_stp.c[j]}: */
                    astp[i].c[j] = (0.80 + 0.20*sin(k+0.3))*adj_stp.c[j];
                  }
              }
          }
      }
    tr2o_debug_params(0, "rad", NI, arad);
    tr2o_debug_params(0, "stp", NI, astp);
  }
   
void tr2o_choose_optimum(uint32_t NI, r2_t pini[], r2_t arad[], r2_t popt[])
  {
    /* Set {popt[0..NI-1]} to {pini[0..NI-1]}, balanced. */ 
    for (uint32_t j = 0;  j < 2; j++)
      { double da[NI];   /* Displacements {popt[i]-pini[i]} along axis {j}. */
        /* Set {da[0..NI-1]} to random displcements; compute their sum: */
        double sumdvar = 0; /* Sum of all {da[0..NI-1]} excluding fixed ones. */
        uint32_t nvar = 0; /* Number of variable coordinates. */
        for (uint32_t i = 0;  i < NI; i++)
          { double ra = arad[i].c[j]; /* Search radius for {p[i]} along axis {j}. */
            if (ra == 0.0)
              { /* Coordinate is fixed: */
                da[i] = 0.0;
              }
            else
              { /* Coordinate is variable: */
                assert(ra > 0.0);
                da[i] = (2*drandom()-1)*ra;
                sumdvar += da[i];
                nvar ++;
              }
          } 
        /* Rebalance {da} to zero sum, add to {pini} to get {popt}: */
        double avg = sumdvar/nvar;
        for (uint32_t i = 0;  i < NI; i++)
          { double ra = arad[i].c[j];
            if (ra > 0.0) { da[i] = da[i] - avg; }
            popt[i].c[j] = pini[i].c[j] + da[i];
          } 
      }
  }

double tr2o_compute_image_mismatch_sqr
  ( uint32_t NI,                       /* Number of images being compared. */
    tr2o_image_eval_proc_t *eval, /* Image evaluator. */
    r2_t p[],       /* Unscaled sampling grid center for each image. */
    i2_t iscale,    /* Image shrink scale in each axis. */
    i2_t wsize,     /* Width of comparison window in each axis (must be odd). */
    double wtx[],   /* Weight table for X displacement. */
    double wty[]    /* Weight table for Y displacement. */
  )  
  { 
    /* Get half-widths of sampling window: */
    demand((wsize.c[0] % 2) == 1, "window width must be odd");
    demand((wsize.c[1] % 2) == 1, "window height must be odd");
    int32_t hwx = (wsize.c[0]-1)/2;
    int32_t hwy = (wsize.c[1]-1)/2;
    
    r2_t scale = (r2_t){{ pow(2.0, iscale.c[0]), pow(2.0, iscale.c[1]) }};
    
    /* Enumerate the sampling points, and get the variances of the values: */
    double sum_wtxy_var = 0;
    double sum_wtxy = 0;
    for (int32_t jy = -hwy; jy <= +hwy; jy++)
      { for (int32_t jx = -hwx; jx <= +hwx; jx++)
          { /* Sample the images at this grid sampling point: */
            double val[NI];
            for (uint32_t i = 0;  i < NI; i++)
              { /* Get samples of image {i}: */
                r2_t *pi = &(p[i]);
                double xsc = pi->c[0]/scale.c[0] + jx;
                double ysc = pi->c[1]/scale.c[1] + jy;
                val[i] = eval(i, iscale, xsc, ysc);
              }
            /* Compute mean and variance, and sampling point weight {wtxy}: */
            double avg, var;
            tr2o_compute_avg_var(NI, val, &avg, &var);
            double wtxy = wtx[jx+hwx]*wty[jy+hwy];
            sum_wtxy_var += wtxy*var;
            sum_wtxy += wtxy;
          }
      }
    return (sum_wtxy == 0 ? 0.0 : sum_wtxy_var / sum_wtxy);
  }

bool_t tr2o_ckeck_result
  ( uint32_t NI, 
    r2_t psol[], 
    double f2sol,
    r2_t pini[], 
    r2_t popt[], 
    r2_t arad[],
    r2_t astp[]
  )
  {
    auto void compare_to_point(char *tag, r2_t pref[]);
      /* Prints the RMS absolute difference {psol[i] - pref[i]}, and 
        the RMS difference {psol[i] - pref[i]} relative to {arad[i]},
        over {i} in {0..NI-1}.  Ignores coordinate {j} of {psol[i],pref[i]} iff
        {arad[i].c[j]} is zero. */
    
    uint32_t indent = 0; /* Indentation */
    fprintf(stderr, "%*schecking the solution\n", indent, "");
    
    bool_t valid = TRUE;
    if (popt != NULL) { compare_to_point("optimum", popt); }
    if (pini != NULL) { compare_to_point("initial", pini); }
    
    fprintf(stderr, "\n");
    return valid;
    
    void compare_to_point(char *tag, r2_t pref[])
      { uint32_t nv = 0;        /* Number of nonzero coordinates in {arad[0..ni-1]}. */
    
        double sum_d2 = 0; /* Total abs squared disp between {psol} and {pref}. */
        double sum_e2 = 0; /* Total rel squared disp between {psol} and {pref}. */

        for (uint32_t i = 0;  i < NI; i++)
          { r2_t *q = &(psol[i]);
            r2_t *o = &(pref[i]);
            r2_t *r = &(arad[i]);
            r2_t d, e;
            r2_sub(q, o, &d);
            for (uint32_t j = 0;  j < 2; j++)
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
        /* Check whether displacement from initial point add to (0,0): */
        fprintf(stderr, "%*s  Abs and rel RMS displacement from %s = %10.6f %10.6f\n", indent, "", tag, rms_d2, rms_e2);
      }
  }

void tr2o_compute_avg_var(uint32_t nz, double z[], double *avgP, double *varP)
  {
    /* Compute the average {avg}: */
    double sum_z = 0; /* Sum of {z[ki]} */
    for (uint32_t i = 0;  i < nz; i++) { sum_z += z[i]; }
    double avg = sum_z/nz; 
    /* Compute the sample variance: */
    double sum_du2 = 0;
    for (uint32_t i = 0;  i < nz; i++) { double dui = z[i] - avg; sum_du2 += dui*dui; }
    double var = sum_du2/nz;
    (*avgP) = avg;
    (*varP) = var;
  }
    
            
    
