/* Tests the cpk_pack routines. */
/* Last edited on 2024-12-31 17:39:03 by stolfi */ 

#define USAGE "test_cpk [-verbose] [-validate] [-plot] [-mag4] INPUT_DIR OUTPUT_DIR"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <affirm.h>
#include <fget.h>
#include <r2.h>
#include <i2.h>
#include <interval.h>

#include <cpk_basic.h>
#include <cpk_valid.h>
#include <cpk_debug.h>
#include <cpk_io.h>
#include <cpk_coords.h>
#include <cpk_main.h>
#include <cpk_main.h>

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);
  
void tcpk_parse_options(int32_t argc, char **argv, cpk_options_t **optP, char **inDirP);
  /* Obtains miscellaneous parameters {*opt} for station allocation,
    either from the command line or by numerical creativity, and 
    the name {inDir} of the directory for input files. */
   
cpk_domain_t *tcpk_get_domain(char *inDir);
  /* Returns a {cpk_domain_t} object that describes the geography
    (including the existing or planned stations). The data is obtained
    from the files 
    
      "{inDir}-urb.txt" -- the urbanized area, as closed polygon(s).
      "{inDir}-mun.txt" -- borders with other municipalities, as polylines. 
      "{inDir}-nat.txt" -- borders with foreign nations, as polylines.
      "{inDir}-exs.txt" -- existing/planned stations, as points.
      "{inDir}-auc.txt" -- predetermined auction points, as points.
      "{inDir}-dis.txt" -- the declared demand points, as points.
    
    Coordinates in the files should be Lon/Lat pairs (in fractional
    degrees). */

cpk_policy_t *tcpk_get_cpk_policy(char *inDir);
  /* Returns a {cpk_policy_t} object that describes the allocation
    policy, including the nominal station range, the radius of auction
    areas, etc.. See the description of {cpk_policy_t}. */

void tcpk_adjust_policy(cpk_policy_t *P, cpk_options_t *opt, bool_t for_demand);
  /* Sets {P->for_demand} as specified, and adjust other parameters (such as
   weight term coefficients) as appropriate. */

r2_vec_t tcpk_all_grid_points(r2_t org, double sx, double sy, uint32_t nx, uint32_t ny);
  /* Generates all points of an hexagonal grid with lower
    left pixel at {org}, {nx} points with step {sx} within each row,
    {ny} rows with step {sy} between rows. Odd-numbered rows
    are offset by {sx/2} in the X direction. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {   
    /* Get command line parameters. */
    /* All file name arguments should be without the extension. */
    /* Input files may be "-" to mean {stdin}. */

    /* Get/fake the miscellaneous options: */
    cpk_options_t *opt;
    char *inDir; /* Directory for input files. */
    tcpk_parse_options(argc, argv, &opt, &inDir);
    
    /* Read the geography data: */
    cpk_domain_t *C = tcpk_get_domain(inDir);

    /* Get/fake the allocation policy parameters: */
    cpk_policy_t *P = tcpk_get_cpk_policy(inDir);
    fprintf(stderr, "dMin = " XY_FFMT "\n", P->dMin);
    fprintf(stderr, "rAuc = " XY_FFMT "\n", P->rAuc);
    fprintf(stderr, "dMun = " XY_FFMT "\n", P->dMun);
    fprintf(stderr, "dNat = " XY_FFMT "\n", P->dNat);

    /* Seed for random numbers: */
    r2_vec_t VJ;
    double WJ;

    /* Initial chronometer reading: */
    double start;
    
    /* Verify clock function: */
    fprintf(stderr, "cpk_cpu_time_1() = %f\n", cpk_cpu_time_1());
    fprintf(stderr, "cpk_cpu_time_2() = %lf\n", cpk_cpu_time_2());
    
    /* Run test with demand: */
    fprintf(stderr, "=== test with demand =========================\n");
    tcpk_adjust_policy(P, opt, TRUE);
    start = cpk_cpu_time_2();
    cpk_choose_auction_points(C, P, opt, &VJ, &WJ);
    fprintf(stderr, "time = %.2f sec\n", (cpk_cpu_time_2() - start)/1.0e6);
    fprintf(stderr, "value = " WT_FFMT "\n", WJ);

    /* Verify clock function: */
    fprintf(stderr, "cpk_cpu_time_1() = %f\n", cpk_cpu_time_1());
    fprintf(stderr, "cpk_cpu_time_2() = %lf\n", cpk_cpu_time_2());
    
    /* Run test without demand: */
    fprintf(stderr, "=== test without demand =========================\n");
    tcpk_adjust_policy(P, opt, FALSE);
    start = cpk_cpu_time_2();
    cpk_choose_auction_points(C, P, opt, &VJ, &WJ);
    fprintf(stderr, "time = %.2f sec\n", (cpk_cpu_time_2() - start)/1.0e6);
    fprintf(stderr, "value = " WT_FFMT "\n", WJ);

    /* Verify clock function: */
    fprintf(stderr, "cpk_cpu_time_1() = %f\n", cpk_cpu_time_1());
    fprintf(stderr, "cpk_cpu_time_2() = %lf\n", cpk_cpu_time_2());
    
    return 0;
  }
  
cpk_domain_t *tcpk_get_domain(char *inDir)
  {
    cpk_domain_t *C;
    C = (cpk_domain_t *)notnull(malloc(sizeof(cpk_domain_t)), "no mem");
    /* Read main polygon: */
    C->Urb = cpk_read_polys(inDir, "urb", TRUE);
    /* Read existing/planned stations: */
    C->Exs = cpk_read_points(inDir, "exs");
    /* Read predetermined auction centers: */
    C->Auc = cpk_read_points(inDir, "auc");
    /* Read municipal and international borders: */
    C->Mun = cpk_read_polys(inDir, "mun", FALSE);
    C->Nat = cpk_read_polys(inDir, "nat", FALSE);
    /* Read declared demand points (DIs): */
    C->Dem = cpk_read_points(inDir, "dis");
    return C;
  }
  
cpk_policy_t *tcpk_get_cpk_policy(char *inDir)
  {
    cpk_policy_t *P;
    P = (cpk_policy_t *)notnull(malloc(sizeof(cpk_policy_t)), "no mem");
    
    /* "Official" application parameters: */
    P->dMin = 4000;     /* Minimum distance between stations. */
    P->rAuc = 1000;     /* Radius of auctioned area. */
    P->dMun = P->rAuc;  /* Min dist of auction ctr to other municipalities. */
    P->dNat = 5000;     /* Min dist of auction ctr to other countries. */

    /* These will be changed later: */
    P->for_demand = FALSE;
    { cpk_wcoeffs_t wc;
      wc.pNUM = 3; wc.cNUM = 1.0;
      wc.pURB = 0; wc.cURB = 1.0;
      wc.pDEM = 1; wc.cDEM = 1.0;
      wc.pCLS = 2; wc.cCLS = 1.0;
      P->wcoeffs = wc;
    }
    return P;
  }
  
void tcpk_parse_options(int32_t argc, char **argv, cpk_options_t **optP, char **inDirP)
  { 
    cpk_options_t *opt;
    opt = (cpk_options_t *)notnull(malloc(sizeof(cpk_options_t)), "no mem");

    /* Default values: */
    char *inDir = ".";
    opt->verbose = FALSE;     /* TRUE prints various diagnostics. */
    opt->validate = FALSE;    /* TRUE runs (expensive) validation checks on solution. */
    opt->plot = FALSE;        /* TRUE plots the solutions. */
    opt->maxSecsGRASP = 30.0; /* Maximum seconds for GRASP optimization. */
    opt->magnify = 1.0;       /* Domain enlargement factor (for testing). */
    opt->seed = 46150317;     /* My favorite random number. */
    opt->outDir = ".";        /* Prefix for output file names. */
    
    /* Parse command line: */
    int32_t i = 1;
    while ((i < argc) && (argv[i][0] == '-'))
      { char *op = argv[i];
        fprintf(stderr, "argv[%d] = %s\n", i, op);
        if (strcmp(op, "-verbose") == 0)
          { opt->verbose = TRUE; }
        else if (strcmp(op, "-validate") == 0)
          { opt->validate = TRUE; }
        else if (strcmp(op, "-plot") == 0)
          { opt->plot = TRUE; }
        else if (strcmp(op, "-mag4") == 0)
          { opt->magnify = 4.0; }
        else 
          { fprintf(stderr, "unrecognized option \"%s\"\n%s\n", op, USAGE); exit(1); }
          i++;
      }
    if (i != argc-2) { fprintf(stderr, "excess/missing arguments\n%s\n", USAGE); exit(1); }
    
    inDir = argv[i];
    opt->outDir = argv[i+1];

    (*optP) = opt;
    (*inDirP) = inDir;
  }

void tcpk_adjust_policy(cpk_policy_t *P, cpk_options_t *opt, bool_t for_demand)
  {
    P->for_demand = for_demand;
    cpk_wcoeffs_t *wc = &(P->wcoeffs);
    if (for_demand)
      { /* Maximize number of auction points, break ties by demand served, etc. */
        wc->pNUM = 0; wc->cNUM = 1.0;
        wc->pDEM = 1; wc->cDEM = 1.0;
        wc->pCLS = 2; wc->cCLS = 1.0;
        wc->pURB = 3; wc->cURB = 1.0;
      }
    else
      { /* Maximize number of auction points, break ties by urban area coverage, etc. */
        wc->pNUM = 0; wc->cNUM = 1.0;
        wc->pURB = 1; wc->cURB = 1.0;
        wc->pDEM = 2; wc->cDEM = 1.0;
        wc->pCLS = 3; wc->cCLS = 1.0;
      }
  }

r2_vec_t tcpk_all_grid_points(r2_t org, double sx, double sy, uint32_t nx, uint32_t ny)
  {
    /* Create the point list: */
    r2_vec_t V = r2_vec_new(nx*ny);

    /* Fill it with a simple hexagonal grid of points: */
    { for (uint32_t iy = 0; iy < ny; iy++)
        { double y = Y(org) + iy*sy;
          for (uint32_t ix = 0; ix < nx; ix++)
            { double x = X(org) + (ix + 0.5*(iy & 1))*sx; 
              V.e[iy*nx + ix] = (r2_t){{x,y}}; 
            }
        }
    }
    return V;
  }
