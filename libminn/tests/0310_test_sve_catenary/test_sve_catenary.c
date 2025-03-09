/* test_sve_catenary --- test of {sve_minn.h} for hanging-chain energy minimization. */
/* Last edited on 2025-02-16 20:40:18 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <values.h>

#include <bool.h>
#include <sign.h>
#include <argparser.h>
#include <jsrandom.h>
#include <jsstring.h>
#include <jsprintf.h>
#include <affirm.h>
#include <rmxn.h>
#include <rn.h>
#include <jsfile.h>
#include <vec.h>

#include <sve_minn.h>
#include <minn_plot.h>

/* GENERAL PARAMETERS */

typedef struct options_t
  { uint32_t nk; /* Number of links in the chain. */
    double wd; /* Ideal distance between attachment points. */
  } options_t;

#define MAXLINKS 256
  /* Max number of links in chain. */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

options_t *get_options(int32_t argc, char **argv);
  /* Parses the command-line options. */

void find_chain_shape(uint32_t nk, double wd);
  /* Tries to find node positions that minimize the 
    potential energy of an {nk}-link chain with unit-length links,
    whose left endpoint is attached to {(0,0)} and whose
    right endpoint slides on a horizontal rail at {Y=0}
    attached by a spring to {(wd,0)}. */ 

void write_solution
  ( char *prefix, 
    char *tag, 
    sve_goal_t *F, 
    uint32_t nc,
    double c[],
    double Fc,
    uint32_t nk,
    double x[],
    double y[],
    double wd
  );
  /* Writes the node positions to stderr and, if {prefix} is not NULL
    also to file "{prefix}-{tag}.dat". Also prints the function value {Fx}
    to {stderr}. */
  
void write_node_positions(char *prefix, char*tag, uint32_t nk, double x[], double y[]);
  /* Writes the chain node positions to a file named "{fname}". Each
    line contains the node index and two values {(x[i],y[i])} for {i =
    0..nk}.  Note that {x} and {y} must have {nk+1} elements each.
    If {fname} is NULL, writes to {stderr} instead. */

void compute_node_positions(uint32_t nk, double c[], double x[], double y[]);
  /* Computes the node positions implied by the coefficients
    {c[0..nk-2]}. Coefficient {c[i]} encodes the angle between links
    {i} and {i+1} of the chain, in the range {{-1 _ +1]}; where {-1}
    means a bend downward by 175 degrees and {+1} means a bend upwards 
    by 175 degrees..
    
    The arrays {x,y} must have {nk+1} elements each. Each pair
    {(x[j],y[j])} will be the coordinates of articulation {j} from the
    left. In particular, the point {(x[0],y[0])} is always {(0,0)} and
    point {(x[nk-1],y[nk-1])} is {(a,0)} for some {a}. */

double chain_energy(uint32_t nk, double x[], double y[], double wd);
  /* The potential energy of the chain, plus the elastic energy of a spring
    that tries to keep the right end-point of the chain at abscissa {x[nk] = 1}. */

void plot_energy
  ( char *prefix,
    char *tag,
    sve_goal_t *F,
    uint32_t nc,
    double c[]
  );
  /* Writes a FNI file "{prefix}-{tag}-plt.fni"
    with a plot of the energy on an arbitrary 2D
    subspace through the point {c}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { options_t *o = get_options(argc, argv);
    srand(4615);  srandom(4615);
    find_chain_shape(o->nk, o->wd);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }
  
void find_chain_shape(uint32_t nk, double wd)
  { 
    /* Output file names: */
    char *prefix = jsprintf("out/%03d", nk);
    
    /* Number of minimization variables: */
    uint32_t nc = nk - 1;
    
    /* Working storage for the goal function: */
    double x[nk+1], y[nk+1];  /* Node positions. */
    
    auto double sve_goal(uint32_t n, double c[]); 
      /* The goal function for optimization. */
      
    double cprev[nc]; /* Guess in previous call of {sve_OK} function. */
    uint32_t nok = 0;      /* Counts iterations (actually, calls to {sve_OK}). */
    bool_t sve_debug = FALSE;
    bool_t sve_debug_probes = FALSE;
    
    auto bool_t sve_OK(uint32_t iter, uint32_t n, double c[], double Fc, double dist, double step, double radius); 
      /* Acceptance criterion function. */

    double c[nc];     /* Initial guess and final solution. */
    
    /* Initialize {c} with a semicircle: */
    for (uint32_t i = 0;  i < nc; i++) { c[i] = (175.0/180.0)/nc; }
    
    /* Print and write the initial guess: */
    fprintf(stderr, "initial guess:\n");
    compute_node_positions(nk, c, x, y);
    double Fc = sve_goal(nc, c);
    write_solution(prefix, "ini", sve_goal, nc, c, Fc, nk, x, y, wd);
    
    /* Optimize iteratively: */
    double *ctr = NULL;
    double dMax = (175.0/180.0)*M_PI;
    bool_t dBox = FALSE;
    double rMin = 0.000001;
    double rMax = 0.50*M_PI;
    double rIni = 0.10*M_PI;
    double minStep = 0.01*rMin;
    uint32_t maxIters = 200;
    sign_t dir = -1;
    
    sve_minn_iterate
      ( nc, &sve_goal, &sve_OK, NULL,
        c, &Fc, dir, 
        ctr, dMax, dBox, rIni, rMin, rMax, 
        minStep, maxIters, sve_debug, sve_debug_probes
      );
    
    /* Print and write final solution: */
    fprintf(stderr, "final solution:\n");
    compute_node_positions(nk, c, x, y);
    write_solution(prefix, "fin", sve_goal, nc, c, Fc, nk, x, y, wd);
    return;
      
    double sve_goal(uint32_t n, double c[])
      { assert(n == nc);
        /* Compute the potential energy of the chain: */
        compute_node_positions(nk, c, x, y);
        double Fc = chain_energy(nk, x, y, wd);
        return Fc;
      }
      
    bool_t sve_OK(uint32_t iter, uint32_t n, double c[], double Fc, double dist, double step, double radius)
      { assert(n == nc);
        fprintf(stderr, "iteration %d\n", nok);
        if (nok > 0)
          { double d = rn_dist(n, cprev, c);
            fprintf(stderr, "change = %16.10f\n", d);
          }
        if (sve_debug) 
          { compute_node_positions(nk, c, x, y);
            write_solution(NULL, "tmp", sve_goal, nc, c, Fc, nk, x, y, wd);
          }
        /* Save guess in {xprev} for next call: */
        rn_copy(n, c, cprev); 
        nok++;
        return FALSE;
      }
    
  }
  
void write_solution
  ( char *prefix, 
    char *tag, 
    sve_goal_t *sve_goal, 
    uint32_t nc,
    double c[],
    double Fc,
    uint32_t nk,
    double x[],
    double y[],
    double wd
  )
  { fprintf(stderr, "  positions =\n");
    write_node_positions(NULL, tag, nk, x, y); 
    fprintf(stderr, "\n");
    
    double E = chain_energy(nk, x, y, wd);
    fprintf(stderr, "  potential = %+24.16e\n", E);
    fprintf(stderr, "\n");
    
    double FcN = sve_goal(nc, c);
    fprintf(stderr, "goal function = %+24.16e %+24.16e\n", Fc, FcN);
    demand(Fc == FcN, "inconsistent function value on return");

    if (prefix != NULL) 
      { write_node_positions(prefix, tag, nk, x, y);
        plot_energy(prefix, tag, sve_goal, nc, c);
      }
  }

void write_node_positions(char *prefix, char *tag, uint32_t nk, double x[], double y[])
  { FILE *wr;
    if (prefix == NULL)
      { wr = stderr; }
    else
      { char *fname = jsprintf("%s-%s.dat", prefix, tag);
        wr = open_write(fname, TRUE);
        free(fname);
      }
    for (uint32_t i = 0;  i <= nk ; i++)
      { fprintf(wr, "%5d %12.8f %12.8f\n", i, x[i], y[i]); }
    if ((wr != stderr) && (wr != stdout)) { fclose(wr); }
  }

void plot_energy
  ( char *prefix,
    char *tag,
    sve_goal_t *sve_goal,
    uint32_t nc,
    double c[]
  )
  { /* Choose two orthogonal deformation modes {ua,ub} with mags {ra,rb}: */
    double va[nc], ua[nc];
    for (uint32_t i = 0;  i < nc; i++) { va[i] = 0.05; }
    double ra = rn_dir(nc, va, ua);
    double vb[nc], ub[nc];
    for (uint32_t i = 0;  i < nc; i++) { vb[i] = i*0.05/nc; }
    double dba = rn_dot(nc, vb, ua);
    rn_mix_in(nc, -dba, ua, vb);
    double rb = rn_dir(nc, vb, ub);
    /* Plot the energy as image: */
    double step = fmax(ra, rb)/30;
    float_image_t *img = minn_plot_2D_float_image(nc, c, ua, ra, ub, rb, TRUE, step, sve_goal);
    char *fname = jsprintf("%s-%s-plt.fni", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, img);
    fclose(wr);
    free(fname);
  }

void compute_node_positions(uint32_t nk, double c[], double x[], double y[])
  { uint32_t nc = nk - 1;
    /* Initial node is always {(0,0)}: */
    x[0] = y[0] = 0;
    /* Assume that link 0 is horizontal: */
    double aabs = 0;        /* Angle of next link, relative to X axis. */
    x[1] = 1; y[1] = 0; 
    /* Compute {(x[i],y[i]) for {i} in {2..nk}: */
    for (uint32_t i = 0;  i < nc; i++)
      { double arel = c[i]*175.0/180.0*M_PI; /* Angle of link relative to previous link. */
        aabs = aabs + arel;
        double dx = cos(aabs), dy = sin(aabs);
        x[i+2] = x[i+1] + dx;
        y[i+2] = y[i+1] + dy;
      }
    /* Compute the rotation needed to bring the last node to the {x}-axis: */
    double r = - atan2(y[nk], x[nk]);
    /* Apply the rotation to all joints: */
    double cr = cos(r), sr = sin(r);
    for (uint32_t i = 1;  i <= nk; i++)
      { double xr = cr*x[i] - sr*y[i];
        double yr = sr*x[i] + cr*y[i];
        x[i] = xr; y[i] = yr;
      }
  }

double chain_energy(uint32_t nk, double x[], double y[], double wd)
  { /* Compute the gravitational energy {G}: */
    double G = 0;
    for (uint32_t i = 0;  i < nk; i++) 
      { double yb = (y[i] + y[i+1])/2;
        G += yb;
      }
    /* Compute elastic energy {E} of spring at right end: */
    double dx = x[nk] - wd;
    double E = dx*dx;
    /* Compute the total energy: */
    return G + E;
  }
  
options_t *get_options(int32_t argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    o->nk = (uint32_t)argparser_get_next_int(pp, 1, MAXLINKS);
    o->wd = argparser_get_next_double(pp, 0.0, DBL_MAX);
    argparser_finish(pp);
    return o;
  }
