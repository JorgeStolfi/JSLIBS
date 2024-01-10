#define PROG_NAME "ppvip"
#define PROG_DESC "the PPV Image Processor for manipulation of PPV images"
#define PROG_VERS "1.0"

/* Copyright © 2006 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2007-10-28 19:50:36 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [-seed NUM] \\\n" \
  "    [VAR=]EXPR \\\n" \
  "    ... \n" \
  "  where EXPR is one of:\n" \
  "    NUM\n" \
  "    VAR\n" \
  "    read([FMT:]FNAME)\n" \
  "    write(EXPR,[FMT:]FNAME)\n" \
  "    add(EXPR,EXPR)\n" \
  "    crop(EXPR,AXIS,NUM,NUM)\n" \
  "  VAR is a variable identifier,\n" \
  "  NUM is a numeric constant,\n" \
  "  AXIS is a nonnegative integer,\n" \
  "  FMT is a file format tag (PPM,PPV,...),\n" \
  "  FNAME is a file name, possibly with escaped chars."

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  ???\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in mar/2006 by Jorge Stolfi, IC-UNICAMP."

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
  
#include <ppv_array.h>
#include <ppv_io.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <math.h>

typedef struct options_t 
  { int seed;            /* Seed for random number generators. */
    command_t *cmds;     /* List of commands. */
  } options_t;

int main (int argc, char **argv);
options_t *get_options(int argc, char **argv);
double parse_double(int *argn, int argc, char **argv, double lo, double hi);
void data_error(int line, char *msg);
void arg_error(char *msg, int argn, char *arg);

int main (int argc, char **argv)
  { 
    options_t *o = get_options(argc, argv);
    
    double c[o->nv*o->nv];

    srandom(o->seed);
    gen_arc_cost_matrix(o->nv, o->distr, o->d, c);

    /* Enumerate tours: */
    vtx_t *p;
    int np;
    enum_perms(o->nv, o->tour, o->npmax, &np, &p);

    /* Compute their costs: */
    double Z[np];
    compute_perm_costs(o->nv, o->tour, c, np, p, Z);
    
    /* Index-sort of costs: */
    int ix[np];
    sort_costs(np, ix, Z);

    /* Fit global sphere-slice model: */
    process_global_model(o, np, ix, p, Z);
    
    /* Fit asymptotic power model: */
    int ns = 2*(int)sqrt(np);
    if (ns > np/2) { ns = np/2; }
    if (ns > o->nsmax) { ns = o->nsmax; }
    process_asymptotic_model(o, ns, +1, np, ix, p, Z);
    process_asymptotic_model(o, ns, -1, np, ix, p, Z);
    
    return 0;
  }

void sort_costs(int np, int *ix, double *Z)
  { 
    int j; 
    for(j = 0; j < np; j++) { ix[j] = j; }
    isrt_heapsort(ix, np, cmp_int, +1);
  }

void process_global_model(options_t *o, int np, int *ix, vtx_t *p, double *Z)
  { 
    /* Fit global model: */
    int Sdim; double Zmid, Zrad;
    int useSdim = o->sdim;
    fit_global_model(np, ix, Z, o->nv, o->tour, useSdim, &Sdim, &Zmid, &Zrad);
    /* Compute ranks {TM} predicted by global model: */
    double TM[np];
    compute_global_model_ranks(np, Z, TM, Sdim, Zmid, Zrad);
    /* Open plot file: */
    FILE *wr = open_plot_file("glob", o, "Global model");
    /* Write model parameters: */
    fprintf(wr, "Sdim = %4d\n", Sdim);
    fprintf(wr, "Zmid = %24.16e\n", Zmid);
    fprintf(wr, "Zrad = %24.16e\n", Zrad);
    /* Write plot data: */
    write_plot_data(wr, o->nv, np, ix, Z, TM, NULL, p);
    fclose(wr);
  }

void process_asymptotic_model(options_t *o, int ns, sign_t dir, int np, int *ix, vtx_t *p, double *Z)
  { 
    char *hilo = (dir > 0 ? "lo" : "hi");
    char *title = NULL;
    asprintf(&title, "asymptotic model for %d out of %d samples at %s end", ns, np, hilo);
    
    fprintf(stderr, "%s: Fitting %s\n", __FUNCTION__, title); 
    int nv = o->nv;
    demand(ns <= np, "bad ns");
    /* Extract subset to analyze: */
    double ZS[ns];
    vtx_t ps[ns*nv];
    int j;
    for (j = 0; j < ns; j++)
      { int k = (dir > 0 ? ix[j] : ix[np - 1 - j]);
        ZS[j] = (dir > 0 ? Z[k] : -Z[k]);
        vtx_t *v = &(p[k*nv]);
        vtx_t *vs = &(ps[j*nv]);
        int i;
        for (i = 0; i < nv; i++) { vs[i] = v[i]; }
      }
    /* Fit asymptotic model: */
    double Expt, Zref, Coef;
    int useExpt = (o->sdim <= 0 ? 0.0 : ((double)o->sdim + 1)/2);
    fit_asymp_model(ns, ZS, np, nv, useExpt, &Expt, &Zref, &Coef);
    /* Compute costs according to asymptotic model: */
    double ZT[ns];
    compute_asymp_model_costs(ns, ZT, np, Expt, Zref, Coef);
    /* Write plot file: */
    char *tag = txtcat("asym-", hilo);
    FILE *wr = open_plot_file(tag, o, title);
    free(tag); free(title);
    /* Write model parameters: */
    fprintf(wr, "Expt = %6.1f\n", Expt);
    fprintf(wr, "Zref = %24.16e\n", Zref);
    fprintf(wr, "Coef = %24.16e\n", Coef);
    /* Write plot data: */
    write_plot_data(wr, o->nv, ns, NULL, ZS, NULL, ZT, ps);
    fclose(wr);
  }

int compute_num_tours(int nv, bool_t sym)
  {
    if (nv == 0)
      { return 0; }
    else if (nv == 1)
      { return 1; }
    else if (nv == 2)
      { return 1; }
    else if (nv == 3)
      { return (sym ? 1 : 2); }
    else
      { int nt1 = compute_num_tours(nv-1, sym);
        if(nt1 > NP_MAX/(nv-1)) { return NP_MAX+1; }
        return nt1*(nv-1); 
      }
  }

int compute_num_spins(int nv)
  {
    if (nv == 0)
      { return 1; }
    else if (nv == 1)
      { return 0; }
    else 
      { int nq1 = compute_num_quasi_spins(nv-1);
        if (nq1 > NP_MAX/(nv-1)) { return NP_MAX+1; }
        return nq1*(nv-1); 
      }
  }

int compute_num_quasi_spins(int nv)
  /* Number of perms of {0..nv-1} which have at least one 
    cycle and no loops except at vertex {0}. */
  {
    if (nv == 0)
      { return 0; }
    else if (nv == 1)
      { return 1; }
    else
      { /* Count quasi-spins that HAVE a fixed point at 0: */
        int ns1 = compute_num_spins(nv-1);
        /* Count quasi-spins that DO NOT HAVE a fixed point at 0: */
        int nq1 = compute_num_quasi_spins(nv-1);
        if (nq1 > (NP_MAX - ns1)/(nv-1)) { return NP_MAX+1; }
        return ns1 + nq1*(nv-1); 
      }
  }

void enum_perms(int nv, bool_t tour, int npmax, int *np, vtx_t **p)
  {
    demand(nv >= 0, "invalid number of vertices");
    demand(nv <= NV_MAX, "too many vertices");
    
    /* Eventually this will be a parameter: */
    bool_t sym = TRUE; /* If {TRUE}, ignores tours that differ by reversal. */

    /* Compute number {np} of perms to generate. */
    if (npmax != 0) 
      { demand(npmax <= NP_MAX, "asking for too many perms");
        (*np) = npmax; 
      }
    else 
      { if (tour)
          { (*np) = compute_num_tours(nv, sym); }
        else
          { (*np) = compute_num_spins(nv); }
        demand((*np) <= NP_MAX, "there are too many perms");
      }
    
    /* Allocate perm table: */
    (*p) = (vtx_t *)notnull(malloc((*np)*nv*sizeof(vtx_t)), "no mem");
    
    vtx_t v[nv];
    
    auto void do_enum(int k);
      /* If {npmax == 0}, generates all valid perms of {v} that start with
        {v[0..k-1]}, else generates a random one of those perms. Either
        way, stores the enumeradted perm(s) in the {p} array starting at
        {p[ngen*nv]}, and increments {ngen}. Upon exit, {v} is
        restored to its original contents. */

    int ngen = 0; /* Counts output perms. */

    void do_enum(int k)
      { 
        if ((!tour) && (k > 0) && (v[k-1] == k-1))
          { /* Skip spins with loops: */
            return;
          }
        else if (k == nv-1) 
          { /* Perm is complete. */
            /* If costs are symmetric, ignore "decreasing" tours: */
            if (tour && sym && (k > 0) && (v[k] < v[1])) { return; }
            /* Store perm in {p} array: */
            assert (ngen < (*np));
            { int j; 
              vtx_t *pp = &((*p)[ngen*nv]); 
              for (j = 0; j < nv; j++) { pp[j] = v[j]; }
            }
            ngen ++; 
          }
        else if (npmax == 0)
          { /* Try all remaining choices for {v[k]}, recurse: */
            int i;
            for (i = k; i < nv; i++)
              { int t = v[k]; v[k] = v[i]; v[i] = t;
                do_enum(k+1);
                t = v[k]; v[k] = v[i]; v[i] = t;
              }
          }
        else
          { /* Pick a random choice for {v[k]}, recurse: */
            int i = k + random() % (nv - k);
            int t = v[k]; v[k] = v[i]; v[i] = t;
            do_enum(k+1);
            t = v[k]; v[k] = v[i]; v[i] = t;
          }
      }

    /* Initial perm (identity): */
    int i;
    for (i = 0; i < nv; i++) { v[i] = i; }

    /* Enumerate desired perms: */
    while (ngen < (*np))
      { /* generate all perms or or one perm: */
        if (tour) 
          { /* Enumerate all tours with {v[0] = 0}: */
            do_enum(1);
          }
        else
          { /* Enumerate all derangements: */
            do_enum(0);
          }
      }

    fprintf(stderr, "generated %d valid perms\n", ngen);
    assert(ngen == (*np));
  }

void write_plot_data(FILE *wr, int nv, int np, int *ix, double *Z, double *TM, double *ZM, vtx_t *p)
  { int j;
    for (j = 0; j < np; j++)
      { /* Index of element with actual rank j: */
        int k = (ix == NULL ? j : ix[j]);
        /* Actual relative rank: */
        double Tk = ((double)j + 0.5)/((double) np);
        /* Actual cost: */
        double Zk = Z[k];
        /* Model-predicted relative rank: */
        double TMk = (TM == NULL ? Tk : TM[k]);
        /* Model-predicted cost: */
        double ZMk = (ZM == NULL ? Zk : ZM[k]);
        /* Write data: */
        fprintf(wr, "%10d", j);
        fprintf(wr, " %9d", k);
        fprintf(wr, " %18.10e", Tk);
        fprintf(wr, " %18.10e", Zk);
        fprintf(wr, " %18.10e", TMk);
        fprintf(wr, " %18.10e", ZMk);
        if (SHOW_PERM) 
          { /* Write permutation: */
            vtx_t *v = &(p[k*nv]); 
            write_perm(wr, nv, v);
          }
        /* Finish line: */
        fprintf(wr, "\n");
      }
    fflush(wr);
  }

void write_perm(FILE *wr, int nv, vtx_t *v)
  {
    int i;
    for (i = 0; i < nv; i++) { fprintf(wr, " %d", v[i]); } 
  }

FILE *open_plot_file(char *tag, options_t *o, char *title)
  {
    char *fname = NULL;
    affirm(asprintf(&fname, "%s-%s.dat", o->outName, tag) > 0, "no mem");
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    
    fprintf(wr, "# %s\n", title);
    write_params(wr, o->seed, o->nv, o->tour, o->npmax, o->distr, o->d);
    return wr;
  }

options_t *get_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    
    int argn = 1;
    
    /* Preliminary defaults: */
    o->seed = 4615;
    o->nv = -1;
    o->tour = TRUE;
    o->npmax = 0;
    o->nsmax = 200;
    o->distr = DISTR_GEN_USPHE;
    o->d = -1;
    o->sdim = 0;
    o->outName = NULL;
    
    /* Parse command line options: */

    while ((argn < argc) && (argv[argn][0] == '-') && (argv[argn][1] != '\0'))
      { if (strcmp(argv[argn], "-seed") == 0)
          { argn++;
            o->seed = (int)parse_double(&argn, argc, argv, 1.0, ((double)INT_MAX));
          } 
        else if (strcmp(argv[argn], "-nv") == 0)
          { argn++;
            o->nv = (int)parse_double(&argn, argc, argv, 1.0, ((double)NV_MAX));
          }
        else if (strcmp(argv[argn], "-tour") == 0)
          { argn++;
            o->tour = ((int)parse_double(&argn, argc, argv, 0.0, 1.0) != 0);
          }
        else if (strcmp(argv[argn], "-npmax") == 0)
          { argn++;
            o->npmax = (int)parse_double(&argn, argc, argv, 0.0, ((double)NP_MAX));
          }
        else if (strcmp(argv[argn], "-nsmax") == 0)
          { argn++;
            o->nsmax = (int)parse_double(&argn, argc, argv, 0.0, ((double)NP_MAX));
          }
        else if (strcmp(argv[argn], "-distr") == 0)
          { argn++;
            o->distr = distr_from_name(argv[argn]); argn++;
            o->d = (int)parse_double(&argn, argc, argv, 0.0, ((double)INT_MAX));
          }
        else if (strcmp(argv[argn], "-sdim") == 0)
          { argn++;
            o->sdim = (int)parse_double(&argn, argc, argv, 0.0, ((double)NV_MAX));
          }
        else if (strcmp(argv[argn], "-outName") == 0)
          { argn++;
            o->outName = argv[argn]; argn++;
          }
        else 
          { arg_error("invalid option", argn, argv[argn]); exit(1); } 
      }

    if (argn != argc) { arg_error("excess arguments", argn, argv[argn]); exit(1); } 
    
    /* Complete defaults: */
    if (o->nv <= 0) 
      { arg_error("must define \"-nv\"", -1, NULL); }
    if (o->outName == NULL) 
      { arg_error("must define \"-outName\"", -1, NULL); }

    argparser_finish(pp);

    return o;
  }

double parse_double(int *argn, int argc, char **argv, double lo, double hi)
  { char *rest;
    double x;
    if (*argn >= argc) 
      { arg_error("argument value is missing", (*argn)-1, argv[(*argn)-1]); }
    x = strtod(argv[*argn], &rest);
    if (*rest != '\0') 
      { arg_error("invalid number", *argn, argv[*argn]); } 
    if ((x < lo) || (x > hi)) 
      { arg_error("out of range", *argn, argv[*argn]); }
    ++(*argn);
    return x;
  }
  
void data_error(int line, char *msg)
  {
    fprintf(stderr, "%s:%d: **%s\n", "-", line, msg);
    exit(1);
  }

void arg_error(char *msg, int argn, char *arg)
  {
    fprintf(stderr, "%s: **", PROG_NAME);
    if ((argn >= 0) && (arg != NULL))
      { fprintf(stderr, " argv[%d] = %s:", argn, arg); }
    fprintf(stderr, " %s\n", msg);
    fprintf(stderr, "usage: %s \\\n%s", PROG_NAME, PROG_HELP);
    exit(1);
  }     
   
