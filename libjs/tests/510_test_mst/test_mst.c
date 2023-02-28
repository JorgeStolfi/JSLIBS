/* Last edited on 2023-02-18 13:21:28 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <epswr.h>
#include <jsrandom.h>
#include <jsstring.h>
#include <r2.h>

#include <mst.h>

int32_t main(int32_t argc, char **argv);

r2_t *makesites(int32_t N, bool_t normal, bool_t verbose);

void plot_mst (int32_t N, r2_t st[], int32_t P[], double C[], bool_t axes, char *prefix);
  /* Writes the plot file named "out/{prefix}-{page}.eps". */

epswr_figure_t *new_eps_figure
  ( char *prefix,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    int32_t ncap
  );
  /* Opens an EPS file called "out/{prefix}.eps" and sets the 
  client plot window to the rectangle {[xmin_xmax]×[ymin_ymax]} in mm.
  Leaves space at the bottom for {ncap} lines of caption. */

void close_eps_stream(epswr_figure_t *eps, bool_t eps_format);
double max_site_radius(int32_t N, r2_t st[]);
void draw_mst_edges(epswr_figure_t *eps, int32_t N, r2_t st[], int32_t P[], double C[]);
void draw_sites(epswr_figure_t *eps, int32_t N, r2_t st[], int32_t P[], double C[]);

int32_t main(int32_t argc, char **argv)
  { int32_t N = atoi(argv[1]);
    char *distr_name = argv[2];
    char *acost_name = argv[3]; /* Distance function. */
    
    bool_t verbose = (N <= 20);
    
    char *prefix = NULL;
    asprintf(&prefix, "mst_%04d_%s_%s", N, distr_name, acost_name);
    
    bool_t normal = (strcmp(distr_name, "normal") == 0);
    mst_arc_cost_t *acost;
    bool_t plot_axes;
    
    /* Choose the cost functions: */
    auto double edist(int32_t i, int32_t j);         /* Euclidean distance between sites. */
    auto double qdist(int32_t i, int32_t j);         /* Euclidean distance with axes barriers. */
    
    if (strcmp(acost_name, "edist") == 0)
      { acost = &edist; plot_axes = FALSE; }
    else if (strcmp(acost_name, "qdist") == 0)
      { acost = &qdist; plot_axes = TRUE; }
    else
      { affirm(FALSE, "invalid arc cost function"); }
    
    /* Pick the sites: */
    r2_t *st = makesites(N, normal, verbose);

    /* Output vectors for {mst_build_complete}: */
    int32_t *P = notnull(malloc(N*sizeof(int32_t)), "no mem"); /* Predecessor map. */
    double *C = notnull(malloc(N*sizeof(double)), "no mem"); /* Cost map. */
    
    mst_build_complete(N, acost, P, C, verbose);

    plot_mst(N, st, P, C, plot_axes, prefix);
    
    free(prefix);
    free(P);
    free(C);
    free(st);
    return(0);
    
    /* Local proc implementations: */
    
    double edist(int32_t i, int32_t j)
      { return r2_dist(&(st[i]), &(st[j])); }
    
    double qdist(int32_t i, int32_t j)
      { r2_t *si = &(st[i]);
        r2_t *sj = &(st[j]);
        double xp = si->c[0]*sj->c[0];
        double yp = si->c[1]*sj->c[1];
        if ((xp <= 0) || (yp <= 0))
          { return +INF; }
        else
          { return r2_dist(si, sj); }
      }
  }

r2_t *makesites(int32_t N, bool_t normal, bool_t print)
  { 
    r2_t *st = (r2_t *)notnull(malloc(N*sizeof(r2_t)), "no mem");
    
    srandom(4615);

    for (int32_t i = 0; i < N; i++) 
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

void plot_mst (int32_t N, r2_t st[], int32_t P[], double C[], bool_t axes, char *prefix)
  { 
    double r = max_site_radius(N, st);
    
    /* Start a new picture: */
    double wm = 2.4*r;
    double xmin = -wm/2; double xmax = +wm/2;
    double ymin = -wm/2; double ymax = +wm/2;
  
    /* Create Postscript document or EPS figure stream. */
    int32_t ncap = 1;
    epswr_figure_t *eps = new_eps_figure(prefix, xmin,xmax, ymin,ymax, ncap);
    
    if (axes)
      { /* Plot coordinate axes: */
        epswr_set_pen(eps, 0,1,1, 0.20f, 0.0, 0.0);
        epswr_axis(eps, epswr_axis_HOR, 0, xmin, xmax);
        epswr_axis(eps, epswr_axis_VER, 0, ymin, ymax);
      }
    
    /* Plot forest edges, solid: */    
    epswr_set_pen(eps, 0,0,0, 0.20f, 0.0, 0.0);
    draw_mst_edges(eps, N, st, P, C);
    
    /* Plot sites: */
    epswr_set_pen(eps, 0,0,0, 0.20f, 0.0, 0.0);
    draw_sites(eps, N, st, P, C);

    /* Add caption and frame: */
    epswr_set_fill_color(eps, 0,0,0);
    epswr_text(eps, "Minimum Spanning Forest", FALSE, 0.5, TRUE, FALSE);
    epswr_set_pen(eps, 0,0,0, 0.20, 0.0, 0.0);
    epswr_frame(eps);
    /* We are done: */
    epswr_end_figure(eps);
  }

void draw_mst_edges(epswr_figure_t *eps, int32_t N, r2_t st[], int32_t P[], double C[])
  {
    for (int32_t i = 0; i < N; i++)
      { int32_t j = P[i];
        assert((j >= 0) && (j < N));
        if (i != j)
          { r2_t *si = &(st[i]);
            r2_t *sj = &(st[j]);
            epswr_segment(eps, si->c[0], si->c[1], sj->c[0], sj->c[1]);
          }
        else
          { assert(C[i] == +INF); }
      }
  }

void draw_sites(epswr_figure_t *eps, int32_t N, r2_t st[], int32_t P[], double C[])
  {
    double rvile, rroot;
    bool_t labels;
    if (N <= 20) 
      { rvile = 2.00; rroot = 2.50; labels = TRUE; }
    else
      { rvile = 0.50; rroot = 1.00; labels = FALSE; }
    if (labels)
      { epswr_set_label_font(eps, "CourierBold", 7.0); }
    for (int32_t i = 0; i < N; i++) 
      { double x = st[i].c[0];
        double y = st[i].c[1];
        double r;
        if (C[i] == +INF)
          { epswr_set_fill_color(eps, 1.00,0.50,1.00); r = rroot; }
        else
          { epswr_set_fill_color(eps, 1.00,1.00,0.50); r = rvile; }
        epswr_dot(eps, x, y, r,  TRUE, TRUE); 
        if (labels)
          { char *lab = NULL;
            asprintf(&lab, "%d", i);
            epswr_set_fill_color(eps, 0,0,0);
            epswr_label(eps, lab, "0", x, y, 0.0, TRUE, 0.5,0.5, TRUE,FALSE);
            free(lab);
          }
      }
  }
  
double max_site_radius(int32_t N, r2_t st[])
  { double r2max = 0.0;
    for (int32_t i = 0; i < N; i++) 
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
    int32_t ncap
  )
  { 
    /* Select a good figure size: */
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
