/* Last edited on 2024-11-20 21:32:58 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <epswr.h>
#include <jsrandom.h>
#include <jsprintf.h>
#include <jsstring.h>
#include <r2.h>

#include <opf.h>

int32_t main(int32_t argc, char **argv);

r2_t *makesites(uint32_t N, bool_t normal, bool_t verbose);

void plot_opf (uint32_t N, r2_t st[], uint32_t P[], char *prefix);
  /* Writes the plot file named "{prefix}.eps". */

epswr_figure_t *new_eps_figure
  ( char *prefix,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    uint32_t ncap
  );

void close_eps_stream(epswr_figure_t *eps);
double max_site_radius(uint32_t N, r2_t st[]);
void draw_opf_edges(epswr_figure_t *eps, uint32_t N, r2_t st[], uint32_t P[]);
void draw_sites(epswr_figure_t *eps, uint32_t N, r2_t st[], uint32_t P[]);

int32_t main(int32_t argc, char **argv)
  { uint32_t N = (uint32_t)atoi(argv[1]);
    char *distr_name = argv[2];
    char *pcost_name = argv[3];
    
    bool_t verbose = (N <= 20);
    
    opf_arc_cost_t *acost;
    opf_path_cost_t *pcost;
    
    /* Choose the cost functions: */
    auto double Dist(uint32_t i, uint32_t j);         /* Euclidean distance between sites. */
    auto double Fmax(double Ca, double Cp); /* Returns {ca+Cp}. */
    auto double Fsum(double Ca, double Cp); /* Returns {ca+Cp}. */
    
    char *prefix = jsprintf("opf_%04d_%s_%s", N, distr_name, pcost_name);

    acost = &Dist;
    if (strcmp(pcost_name, "fmax") == 0)
      { pcost = &Fmax; }
    else if (strcmp(pcost_name, "fsum") == 0)
      { pcost = &Fsum; }
    else
      { affirm(FALSE, "invalid path cost function"); }
    
    /* Pick the sites: */
    bool_t normal = (strcmp(distr_name, "normal") == 0);
    r2_t *st = makesites(N, normal, verbose);

    /* Output vectors for {opf_build_complete}: */
    uint32_t *P = notnull(malloc(N*sizeof(uint32_t)), "no mem"); /* Predecessor map. */
    uint32_t *R = notnull(malloc(N*sizeof(uint32_t)), "no mem"); /* Root map. */
    double *C = notnull(malloc(N*sizeof(double)), "no mem"); /* Cost map. */
    
    /* Pick three sites as roots: */
    for (uint32_t i = 0;  i < N; i++) { C[i] = +INF; }
    if (N >= 1) { C[0] = 0.0; }
    if (N >= 2) { C[1] = 0.0; }
    if (N >= 3) { C[2] = 0.0; }
    
    opf_build_complete(N, C, acost, pcost, P, R, verbose);

    plot_opf(N, st, P, prefix);
    
    free(prefix);
    free(P);
    free(R);
    free(C);
    free(st);
    return(0);
    
    /* Local proc implementations: */
    
    double Dist(uint32_t i, uint32_t j)
      { return r2_dist(&(st[i]), &(st[j])); }
      
    double Fsum(double Ca, double Cp)
      { return Ca + Cp; }
      
    double Fmax(double Ca, double Cp)
      { return fmax(Ca, Cp); }
  }

r2_t *makesites(uint32_t N, bool_t normal, bool_t print)
  { 
    r2_t *st = (r2_t *)notnull(malloc(N*sizeof(r2_t)), "no mem");
    
    srandom(4615);

    for (uint32_t i = 0;  i < N; i++) 
      { if (normal) 
          { /* Gaussian distribution with total variance 2: */
            st[i].c[0] = dgaussrand();
            st[i].c[1] = dgaussrand();
          }
        else 
          { /* uniform distribution: */
            st[i].c[0] = 2*drandom() - 1;
            st[i].c[1] = 2*drandom() - 1;
          }
        if (print) { fprintf(stderr, "%3d (%12.8f, %12.8f)\n", i, st[i].c[0], st[i].c[1]); }
      }
    return st;
  }

void plot_opf (uint32_t N, r2_t st[], uint32_t P[], char *prefix)
  { 
    double r = max_site_radius(N, st);
    
    /* Start a new picture: */
    double wm = 2.4*r;
    double xmin = -wm/2; double xmax = +wm/2;
    double ymin = -wm/2; double ymax = +wm/2;
  
    /* Create Postscript document or EPS figure stream. */
    uint32_t ncap = 1;
    epswr_figure_t *eps = new_eps_figure(prefix, xmin,xmax, ymin, ymax, ncap);
    
    /* Plot forest edges, solid: */    
    epswr_set_pen(eps, 0,0,0, 0.20f, 0.0, 0.0);
    draw_opf_edges(eps, N, st, P);
    
    /* Plot sites: */
    epswr_set_pen(eps, 0,0,0, 0.20f, 0.0, 0.0);
    draw_sites(eps, N, st, P);

    /* Add caption and frame: */
    epswr_set_pen(eps, 0,0,0, 0.10, 0.0, 0.0);
    epswr_text(eps, "Optimum Path Forest", FALSE, 0.5, TRUE, FALSE);
    epswr_set_pen(eps, 0,0,0, 0.20, 0.0, 0.0);
    epswr_frame(eps);
    /* We are done: */
    epswr_end_figure(eps);
  }

void draw_opf_edges(epswr_figure_t *eps, uint32_t N, r2_t st[], uint32_t P[])
  {
    for (uint32_t i = 0;  i < N; i++)
      { uint32_t j = P[i];
        assert((j >= 0) && (j < N));
        if (i != j)
          { r2_t *si = &(st[i]);
            r2_t *sj = &(st[j]);
            epswr_segment(eps, si->c[0], si->c[1], sj->c[0], sj->c[1]);
          }
      }
  }

void draw_sites(epswr_figure_t *eps, uint32_t N, r2_t st[], uint32_t P[])
  {
    double rvile, rroot;
    bool_t labels;
    if (N <= 20) 
      { rvile = 2.00; rroot = 2.50; labels = TRUE; }
    else
      { rvile = 0.50; rroot = 1.00; labels = FALSE; }
    if (labels)
      { epswr_set_label_font(eps, "CourierBold", 7.0); }
    for (uint32_t i = 0;  i < N; i++) 
      { double x = st[i].c[0];
        double y = st[i].c[1];
        double r;
        if (P[i] == i)
          { epswr_set_fill_color(eps, 1.00,0.50,1.00); r = rroot; }
        else
          { epswr_set_fill_color(eps, 1.00,1.00,0.50); r = rvile; }
        epswr_dot(eps, x, y, r,  TRUE, TRUE); 
        if (labels)
          { char *lab = jsprintf("%d", i);
            epswr_set_fill_color(eps, 0,0,0);
            epswr_label(eps, lab, "0", x, y, 0.0, TRUE, 0.5,0.5, TRUE,FALSE);
            free(lab);
          }
      }
  }
  
double max_site_radius(uint32_t N, r2_t st[])
  { double r2max = 0.0;
    for (uint32_t i = 0;  i < N; i++) 
      { double x = st[i].c[0];
        double y = st[i].c[1];
        double r2i = x*x + y*y;
        if (r2i > r2max) { r2max = r2i; }
      }
    return sqrt(r2max);
  }
  
epswr_figure_t *new_eps_figure
  ( char *prefix,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    uint32_t ncap
  )
  { 
    double fontHeight = 10.0;
    double pt_per_mm = (72.0/25.4); /* One mm in pt. */
    double xfigsz = 150.00*pt_per_mm; /* Figure X size excluding margin (pt). */
    double yfigsz = 150.00*pt_per_mm; /* Figure Y size excluding margin (pt). */
    double mrg = 4.0; /* Figure margin width (pt). */
    double mrg_bot = mrg + (ncap == 0 ? 0 : ncap*fontHeight + mrg);

    char *suffix = NULL;
    bool_t eps_verbose = TRUE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( "out", prefix, NULL, -1, suffix, 
        xfigsz, yfigsz, mrg, mrg, mrg_bot, mrg, 
        eps_verbose
      );
    epswr_set_client_window(eps, xmin,xmax, ymin,ymax);
      
    epswr_set_text_geometry(eps, FALSE, 0,xfigsz,  mrg-mrg_bot,-mrg, 0.0); 
    epswr_set_text_font(eps, "Courier", fontHeight);
    return eps;
  }  
