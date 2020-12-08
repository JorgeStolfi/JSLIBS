/* test_sve_charges --- test of {sve_minn.h} with Rutherford's atom potential.  */
/* Last edited on 2017-03-13 22:20:43 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <values.h>
#include <jsstring.h>

#include <bool.h>
#include <sign.h>
#include <jsrandom.h>
#include <argparser.h>
#include <affirm.h>
#include <r3.h>
#include <r3_extra.h>
#include <rmxn.h>
#include <rmxn_extra.h>
#include <rn.h>
#include <jsfile.h>
#include <js.h>
#include <vec.h>
#include <minn_plot.h>

#include <sve_minn.h>

/* GENERAL PARAMETERS */

typedef struct options_t
  { int nq;   /* Number of electrons. */
    int dim;  /* Free dimensions. */
    int sym;  /* Symmetry class (0 = no symmetry). */
  } options_t;

#define MAXNQ 100
  /* Max number of electrons. */

typedef enum {
    symclass_plain,       /* Independent and unconstrained. */
    symclass_centersym,   /* Centrally symmetric. */
    symclass_regular,     /* Corners of regular polygon/polyhedron. */
    symclass_regcenter,   /* Corners of regular polygon/polyhedron plus central particle. */
    symclass_antiprism,   /* Corners of a regular antiprism. */
    symclass_prism,       /* Corners of a regular prism. */
    symclass_prismcenter, /* Corners of a regular prism. */
    symclass_LIMIT        /* First invalid class. */
  } symclass_t;

/* INTERNAL PROTOTYPES */

int main (int argc, char **argv);

options_t *get_options(int argc, char **argv);
  /* Parses the command-line options. */

double electric_potential(int nq, double Q[], double C, r3_t p[]);
  /* The electric potential of a set of pointlike `electrons' with
    charges {Q[0..nq-1]} at positions {p[0..nq-1]}, immersed in a
    uniform spherical `nuclear cloud' with center at the origin whose
    total charge, up to unit distance from the origin, is {C}.
    
    Assumes that the electric potential inside the cloud grows like
    {r^2} whil the potential between point charges decays like {1/r}. */

void find_electron_positions(int nq, int dim, int sym);
  /* Computes electron coordinates {p[0..nq-1]} that minimize the value of
     {electric_potential(nq, Q, C, p)}. Assumes all point electrons
     have charge {Q[i] = -1} and the cloud charge {C = +nq}.  The 
     electron positions are constained to have the last
     {3-dim} coordinates equal to zero, and symmetry class {sym} */

void write_electron_positions(char *fname, int nq, r3_t p[]);
  /* Writes the electron coordinates {p[0..nq-1]} to a file named
    "{fname}". Each line contains an electron index {i} and three coordinates
    of electron {i}. If {fname} is NULL, writes to {stderr} instead. */

void plot_potential
  ( char *fname, 
    int nq, 
    double Q[], 
    double C, 
    int dim, 
    int sym, 
    int nx, 
    double x0[], 
    double R, 
    int N
  );
  /* Writes a file with sampled values of
    {electric_potential(C,nq,Q,dim,p)} for various sample
    configurations {p}, in a format suitable for the {splot} command
    of {gnuplot}.
    
    Each sample configuration {p[0..nq-1]} is obtained from a sample
    parameter vector {x[0..nx-1]}, where {nx = num_parameters(nq, dim,
    sym)}. If {nx == 1}, the parameter vectors are obtained by setting
    {x[0]} to {2*N+1} equally spaced values spanning the
    interval{x0[0] + [-R_+R]}. If {nx >= 2), the sample vectors {x}
    form a square two-dimensional grid in {R^nx} with {2*N+1} nodes
    and side {2*R} on each axis, centered at the configuration
    {x0[0..nx-1]}. If {nx==2}, the two grid axes are {x[0]} and
    {x[1]}; otherwise they are two random orthogonal directions in
    {R^{nx}}. */

/* PARAMETER MAPPING */

int num_parameters(int nq, int dim, int sym);
  /* Returns the number {nx} of optimization parameters, given the charge count {nq},
    the space dimension {dim}, and the symmetry class {sym}. */

void params_to_coords(int nx, double x[], int nq, int dim, int sym, r3_t p[]);
  /* Maps the parameters {x[0..nx-1]} to the charge coordinates {p[0..nq-1].c[0..dim-1]},
    according to the symmetry class {sym}. */

void coords_to_params(int nq, int dim, int sym, r3_t p[], int nx, double x[]);
  /* Maps the charge coordinates {p[0..nq-1].c[0..dim-1]}, assumed
    to have symmetry class [sym}, to the parameters {x[0..nx-1]}. */

void params_to_coords_plain(int nx, double x[], int nq, int dim, r3_t p[]);
void coords_to_params_plain(int nq, int dim, r3_t p[], int nx, double x[]);
  /* Plain mapping. Parameters are the coordinates. */

void params_to_coords_centersym(int nx, double x[], int nq, int dim, r3_t p[]);       
void coords_to_params_centersym(int nq, int dim, r3_t p[], int nx, double x[]);
  /* Mapping for center-symmetric configurations. Parameters are the coordinates
    of the first {nq/2} charges. */

void params_to_coords_regular(int nx, double x[], int nq, int dim, r3_t p[]);       
void coords_to_params_regular(int nq, int dim, r3_t p[], int nx, double x[]);       
  /* Mapping for regular polygons ({dim==2}) or polyhedra ({dim ==3}). 
    The single parameter is the distance from center to the electrons. */

void params_to_coords_regcenter(int nx, double x[], int nq, int dim, r3_t p[]);         
void coords_to_params_regcenter(int nq, int dim, r3_t p[], int nx, double x[]);     
  /* Mapping for regular polygons ({dim==2}) or polyhedra ({dim ==3})
    with {nq-1} corners, plus one particle at the origin. 
    The single parameter is the distance from center to the electrons. */         

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { options_t *o = get_options(argc, argv);
    find_electron_positions(o->nq, o->dim, o->sym);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }
  
void find_electron_positions(int nq, int dim, int sym)
  { 
    /* Output file names: */
    char *prefix = NULL; /* Prefix for output file names. */
    asprintf(&prefix, "out/test-%03d-%1d-%d", nq, dim, sym);
    char *solfname = txtcat(prefix, "-sol.dat"); /* Name of final solution file. */
    char *potfname = txtcat(prefix, "-pot.dat"); /* Name of potential field file. */
    
    int i; 
    /* Local storage: */
    double Q[nq];     /* Signed charge of each electron. */
    double C;         /* Total signed charge in nuclear cloud. */
    r3_t p[nq];       /* Electron positions. */
    int nx = num_parameters(nq, dim, sym);  /* Size of parameter vector. */
    
    int neval = 0; /* Number of function evaluations: */

    auto double F(int n, double x[]); 
      /* The goal function for uptimization. */
      
    double F(int n, double x[])
      { assert(n == nx);
        params_to_coords(n, x, nq, dim, sym, p);
        double Fx = electric_potential(nq, Q, C, p);
        neval++;
        return Fx;
      }
      
    double Xprev[nx]; /* Guess in previous call of {OK} function. */
    int nok = 0;      /* Counts iterations (actually, calls to {OK}). */
    
    auto bool_t OK(int n, double x[], double Fx); 
      /* Acceptance criterion function. */
      
    bool_t OK(int n, double x[], double Fx)
      { assert(n == nx);
        params_to_coords(n, x, nq, dim, sym, p);
        fprintf(stderr, "iteration %d\n", nok);
        if (nok > 0)
          { double d = rn_dist(n, Xprev, x);
            fprintf(stderr, "displacement = %16.10f\n", d);
          }
        fprintf(stderr, "electron positions =\n");
        write_electron_positions(NULL, nq, p);
        /* Save guess in {Xprev} for next call: */
        rn_copy(n, x, Xprev); 
        nok++;
        fprintf(stderr, "potential = %16.10f\n", Fx);
        fprintf(stderr, "\n");
        return FALSE;
      }
      
    /* Define the electron charges {Q[0..nq-1]} and the total cloud charge {C}: */
    C = 0;
    for (i = 0; i < nq; i++) { Q[i] = -1; C += (-Q[i]); }
    
    double x[nx];     /* Initial guess and final solution. */

    /* Initial electron coordinates: */
    srand(4615);  srandom(4615);
    for (i = 0; i < nq; i++) 
      { p[i] = (r3_t){{ 0,0,0 }};
        rn_throw_ball(dim, p[i].c);
        r3_scale(0.80, &(p[i]), &(p[i]));
      }
    fprintf(stderr, "initial coordinates:\n");
    write_electron_positions(NULL, nq, p);
    coords_to_params(nq, dim, sym, p, nx, x);
    
    /* Optimize iteratively: */
    double dMax = 2.000; /* If all electrons are on the surface. */
    bool_t dBox = TRUE; /* Let charges roam in box. */
    double rMin = 0.050;
    double rMax = 0.500;
    double rIni = 0.125;
    double stop = 0.01*rMin;
    int maxIters = 300;
    sign_t dir = -1;
    bool_t debug = FALSE;
    
    double Fx = F(nx, x);
    sve_minn_iterate(nx, &F, &OK, x, &Fx, dir, dMax, dBox, rIni, rMin, rMax, stop, maxIters, debug);
    fprintf(stderr, "iterations = %d\n", nok);
    fprintf(stderr, "function evaluations = %d\n", neval);
    
    /* Ealuate (and print) the final potential: */
    double FxN = F(nx, x);
    fprintf(stderr, "final potential = %16.10f %16.10f\n", Fx, FxN);
    demand(Fx == FxN, "inconsistent function value on return");
    
    /* Print final electron positions and write file for plotting: */
    params_to_coords(nx, x, nq, dim, sym, p);
    fprintf(stderr, "final coordinates:\n");
    write_electron_positions(NULL, nq, p);
    write_electron_positions(solfname, nq, p);

    /* Dump the potential field around the optimum: */
    double RPlot = 1.5;  /* Half-width of plot domain. */
    int NPlot = 40;
    plot_potential(potfname, nq, Q, C, dim, sym, nx, x, RPlot, NPlot);
    
    free(potfname);
    free(solfname);
    free(prefix);
  }

void plot_potential
  ( char *fname, 
    int nq, 
    double Q[], 
    double C, 
    int dim, 
    int sym, 
    int nx, 
    double x0[], 
    double R, 
    int N
  )
  { demand(nx >= 1, "too few coordinates to plot");
    assert(nx == num_parameters(nq, dim, sym));

    r3_t p[nq];  /* Work area: Electron positions. */

    auto double F(int n, double x[]); 
      /* The goal function for uptimization. */
      
    double F(int n, double x[])
      { assert(n == nx);
        params_to_coords(n, x, nq, dim, sym, p);
        double Fx = electric_potential(nq, Q, C, p);
        return Fx;
      }
    
    int NS = 30;
    FILE *wr = open_write(fname, TRUE);
    minn_plot_2D_gnuplot(wr, nx, F, x0, NS, R);
    if ((wr != stderr) && (wr != stdout)) { fclose(wr); }
  }

void write_electron_positions(char *fname, int nq, r3_t p[])
  { FILE *wr;
    if (fname == NULL) 
      { wr = stderr; }
    else
      { wr = open_write(fname, TRUE);
        fprintf(wr, "# fields: index X Y Z dorg\n");
      }
    int i, j, k;
    for (i = 0; i < nq; i++)
      { fprintf(wr, "%5d", i);
        double d2 = 0;
        for (k = 0; k < 3; k++)
          { double Xik = p[i].c[k];
            fprintf(wr, " %12.8f", Xik);
            d2 += Xik*Xik;
          }
        fprintf(wr, " %12.8f", sqrt(d2));
        fprintf(wr, "\n");
      }

    /* Write electron-electron distances as comments: */
    fprintf(wr, "# fields: i j dij\n");
    for (i = 0; i < nq; i++)
      { for (j = 0; j < i; j++)
          { fprintf(wr, "# %5d %5d", i, j);
            double dij = r3_dist(&(p[i]), &(p[j]));
            fprintf(wr, " %12.8f", dij);
            fprintf(wr, "\n");
          }
      }
    
    fflush(wr);
    if ((wr != stderr) && (wr != stdout)) { fclose(wr); }
  }

double electric_potential(int nq, double Q[], double C, r3_t p[])
  { int i, j;
    double Pcloud = 0; /* Potential electrons-vs-cloud. */
    double Ppairs = 0; /* Potential electrons-vs-electrons. */
    /* Fudge factors to avoid infinities and nans: */
    double relfudge = 1.0e-6;  /* Relative fudge factor. */
    double absfudge = 1.0e-12; /* Absolute fudge term. */
    for (i = 0; i < nq; i++)
      { /* Compute the potential {Pi} due to the cloud on electron {i}: */
        /* {Pi} is 0 at center, {-C*Q[i]} at surface (ri = 1): */
        double Pi = 0.0;
        double ri2 = r3_norm_sqr(&(p[i]));
        if (ri2 <= 1.0)
          { /* The charge below electron {i} is {C*ri^3}, acts as if concentrated at origin: */
            Pi = -C*ri2*Q[i];
          }
        else
          { /* The charge below electron {i} is {C}, acts as if concentrated at origin: */
            double ri = sqrt(ri2);
            Pi = -C*Q[i]*(2 - 1/ri);
          }
        /* fprintf(stderr, "  ri =  %8.4f Pi =  %8.4f\n", r, Pi); */
        Pcloud += Pi;
        /* Add the potential of electrons {0..i-1} on electron {i}: */
        for (j = 0; j < i; j++) 
          { double dij2 = r3_dist_sqr(&(p[i]), &(p[j]));
            double rj2 = r3_norm_sqr(&(p[j]));
            double dij = sqrt(dij2 + relfudge*(ri2+rj2) + absfudge*absfudge);
            double Pij = Q[i]*Q[j]/dij;
            /* fprintf(stderr, "  dij = %8.4f Pij = %8.4f\n", d, Pij); */
            Ppairs += Pij;
          }
      }
    /* fprintf(stderr, "  Pcloud =  %8.4f Ppairs =  %8.4f\n", Pcloud, Ppairs); */
    return Pcloud + Ppairs;
  }
  
options_t *get_options(int argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    o->nq = (int)argparser_get_next_int(pp, 1, MAXNQ);
    o->dim = (int)argparser_get_next_int(pp, 1, 3);
    o->sym = (int)argparser_get_next_int(pp, 0, symclass_LIMIT - 1);
    argparser_finish(pp);
    return o;
  }

int num_parameters(int nq, int dim, int sym)
  { switch(sym)
      { 
      case symclass_plain:
        /* No symmetry: */
        return nq*dim;
      case symclass_centersym:
        /* Central symmetry: */
        return dim*(nq / 2);
        
      case symclass_regular:
        /* Regular polyhedron/polygon: */
        if (dim == 3)
          { demand((nq<=4) || (nq==6) || (nq==8) || (nq==12) || (nq==20), "no such regular polyhedron"); }
        return 1;
        
      case symclass_regcenter:
        /* Regular polyhedron/polygon with central body: */
        if (dim == 3) 
          { demand((nq<=5) || (nq==7) || (nq==9) || (nq==13) || (nq==21), "no such regular polyhedron"); }
        return 1;
        
      case symclass_antiprism:
        /* Antiprism: */
        demand(dim == 3, "invalid dim/sym combination");
        return 2;
        
      case symclass_prism:
        /* Prism: */
        demand(dim == 3, "invalid dim/sym combination");
        return 2;
        
      case symclass_prismcenter:
        /* Prism with central body: */
        demand(dim == 3, "invalid dim/sym combination");
        return 2;
        
      default:
        demand(FALSE, "invalid dim/sym combination");
        return 0;
      }
  }

void params_to_coords(int nx, double x[], int nq, int dim, int sym, r3_t p[])
  { switch(sym)
      { 
      case symclass_plain:
        params_to_coords_plain(nx, x, nq, dim, p); break;

      case symclass_centersym:
        params_to_coords_centersym(nx, x, nq, dim, p); break;
  
      case symclass_regular:
        params_to_coords_regular(nx, x, nq, dim, p); break;
        
      case symclass_regcenter:
        params_to_coords_regcenter(nx, x, nq, dim, p); break;
        
      case symclass_antiprism:
      case symclass_prism:
      case symclass_prismcenter:
        /* Prism with central body: */
        demand(FALSE, "unimplemented dim/sym combination");
        break;
        
      default:
        demand(FALSE, "invalid dim/sym combination");
        break;
      }
  }

void coords_to_params(int nq, int dim, int sym, r3_t p[], int nx, double x[])
  { switch(sym)
      { 
      case symclass_plain:
        coords_to_params_plain(nq, dim, p, nx, x); break;

      case symclass_centersym:
        coords_to_params_centersym(nq, dim, p, nx, x); break;
  
      case symclass_regular:
        coords_to_params_regular(nq, dim, p, nx, x); break;
  
      case symclass_regcenter:
        coords_to_params_regcenter(nq, dim, p, nx, x); break;
  
      case symclass_antiprism:
      case symclass_prism:
      case symclass_prismcenter:
        /* Prism with central body: */
        demand(FALSE, "unimplemented dim/sym combination");
        break;
        
      default:
        demand(FALSE, "invalid dim/sym combination");
        break;
      }
  }


void params_to_coords_plain(int nx, double x[], int nq, int dim, r3_t p[])
  { assert(nx == nq*dim);
    int k = 0;
    int i, j;
    for (i = 0; i < nq; i++)
      { p[i] = (r3_t){{ 0,0,0 }};
        for (j = 0; j < dim; j++)
          { p[i].c[j] = x[k]; k++; }
      }
  }

void coords_to_params_plain(int nq, int dim, r3_t p[], int nx, double x[])
  { assert(nx == nq*dim);
    int k = 0;
    int i, j;
    for (i = 0; i < nq; i++)
      { for (j = 0; j < dim; j++) 
          { x[k] = p[i].c[j]; k++; }
      }
  }

void params_to_coords_centersym(int nx, double x[], int nq, int dim, r3_t p[])        
  { assert(nx == dim*(nq/2));
    int k = 0;
    int i, j;
    for (i = 0; i < nq/2; i++)
      { p[i] = (r3_t){{ 0,0,0 }};
        for (j = 0; j < dim; j++) { p[i].c[j] = x[k]; k++; }
        r3_scale(-1, &(p[i]), &(p[nq-1-i]));
      }
    if ((nq % 2) > 0) { p[nq/2] = (r3_t){{ 0,0,0 }}; }
  }

void coords_to_params_centersym(int nq, int dim, r3_t p[], int nx, double x[])        
  { assert(nx == dim*(nq/2));
    int k = 0;
    int i, j;
    for (i = 0; i < nq/2; i++)
      { for (j = 0; j < dim; j++) 
          { x[k] = p[i].c[j] ; k++; }
      }
  }

void params_to_coords_regular(int nx, double x[], int nq, int dim, r3_t p[])        
  { assert(nx == 1);
    double R = x[0];
    if ((dim == 2) || (nq < 4))
      { /* Regular {nq}-gon on the {Z=0} plane: */
        r3_cylindrical_grid(0.0, R, 1, nq, 0.0, p);
      }
    else
      { if (nq == 4)
          { r3_tetrahedron_vertices(R, nq, p); }
        else if (nq == 6)
          { r3_octahedron_vertices(R, nq, p); }
        else if (nq == 8)
          { r3_hexahedron_vertices(R, nq, p); }
        else if (nq == 12)
          { r3_icosahedron_vertices(R, nq, p); }
        else if (nq == 20)
          { r3_dodecahedron_vertices(R, nq, p); }
        else
          { demand(FALSE, "bad nq/dim/sym combination"); }
      }
  }

void coords_to_params_regular(int nq, int dim, r3_t p[], int nx, double x[])        
  { assert(nx == 1);
    assert(nq >= 1);
    x[0] = r3_norm(&(p[0]));
  }

void params_to_coords_regcenter(int nx, double x[], int nq, int dim, r3_t p[])        
  { assert(nx == 1);
    params_to_coords_regular(nx, x, nq-1, dim, p);
    p[nq-1] = (r3_t){{ 0,0,0 }};
  }

void coords_to_params_regcenter(int nq, int dim, r3_t p[], int nx, double x[])        
  { assert(nx == 1);
    assert(nq >= 1);
    x[0] = r3_norm(&(p[0]));
  }
