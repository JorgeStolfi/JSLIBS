/* Last edited on 2022-10-18 20:59:49 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <pswr.h>
#include <jsrandom.h>
#include <jsstring.h>
#include <r2.h>

#include <mst.h>

int main(int argc, char **argv);

r2_t *makesites(int N, bool_t normal, bool_t verbose);

void plot_mst (int N, r2_t st[], int P[], double C[], bool_t axes, char *prefix, bool_t eps_format);
  /* Writes the plot file named "{prefix}-doc.ps" or "{prefix}-{page}.eps". */

PSStream *open_ps_stream(double r, char *prefix, bool_t eps_format);
void close_ps_stream(PSStream *fps, bool_t eps_format);
double max_site_radius(int N, r2_t st[]);
void draw_mst_edges(PSStream *fps, int N, r2_t st[], int P[], double C[]);
void draw_sites(PSStream *fps, int N, r2_t st[], int P[], double C[]);

int main(int argc, char **argv)
  { int N = atoi(argv[1]);
    bool_t normal = (strcmp(argv[2], "normal") == 0);
    char *dname = argv[3]; /* Distance function. */
    char *prefix = argv[4];
    
    bool_t verbose = (N <= 20);
    
    mst_arc_cost_t *acost;
    bool_t plot_axes;
    
    /* Choose the cost functions: */
    auto double edist(int i, int j);         /* Euclidean distance between sites. */
    auto double qdist(int i, int j);         /* Euclidean distance with axes barriers. */
    
    if (strcmp(dname, "edist") == 0)
      { acost = &edist; plot_axes = FALSE; }
    else if (strcmp(dname, "qdist") == 0)
      { acost = &qdist; plot_axes = TRUE; }
    else
      { affirm(FALSE, "invalid arc cost function"); }
    
    /* Pick the sites: */
    r2_t *st = makesites(N, normal, verbose);

    /* Output vectors for {mst_build_complete}: */
    int *P = notnull(malloc(N*sizeof(int)), "no mem"); /* Predecessor map. */
    double *C = notnull(malloc(N*sizeof(double)), "no mem"); /* Cost map. */
    
    mst_build_complete(N, acost, P, C, verbose);

    bool_t eps_format = TRUE;
    plot_mst(N, st, P, C, plot_axes, prefix, eps_format);
    
    free(P);
    free(C);
    free(st);
    return(0);
    
    /* Local proc implementations: */
    
    double edist(int i, int j)
      { return r2_dist(&(st[i]), &(st[j])); }
    
    double qdist(int i, int j)
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

r2_t *makesites(int N, bool_t normal, bool_t print)
  { 
    r2_t *st = (r2_t *)notnull(malloc(N*sizeof(r2_t)), "no mem");
    
    srandom(4615);

    int i;
    for (i = 0; i < N; i++) 
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

void plot_mst (int N, r2_t st[], int P[], double C[], bool_t axes, char *prefix, bool_t eps_format)
  { 
    double r = max_site_radius(N, st);
  
    /* Create Postscript document or EPS figure stream. */
    PSStream *fps = open_ps_stream(r, prefix, eps_format);
    
    /* Start a new picture: */
    double wm = 2.4*r;
    double xmin = -wm/2; double xmax = +wm/2;
    double ymin = -wm/2; double ymax = +wm/2;
    pswr_new_picture(fps, xmin,xmax, ymin, ymax);
    
    if (axes)
      { /* Plot coordinate axes: */
        float penwd = (eps_format ? 0.20f : 0.10f);
        pswr_set_pen(fps, 0,1,1, penwd, 0.0, 0.0);
        pswr_axis(fps, HOR, 0, xmin, xmax);
        pswr_axis(fps, VER, 0, ymin, ymax);
      }
    
    /* Plot forest edges, solid: */    
    float penwd = (eps_format ? 0.20f : 0.10f);
    pswr_set_pen(fps, 0,0,0, penwd, 0.0, 0.0);
    draw_mst_edges(fps, N, st, P, C);
    
    /* Plot sites: */
    pswr_set_pen(fps, 0,0,0, penwd, 0.0, 0.0);
    draw_sites(fps, N, st, P, C);

    if (! eps_format)
      { /* Add caption and frame: */
        pswr_set_pen(fps, 0,0,0, 0.10, 0.0, 0.0);
        pswr_add_caption(fps, "Minimum Spanning Forest", 0.5);
        pswr_set_pen(fps, 0,0,0, 0.20, 0.0, 0.0);
        pswr_frame(fps);
      }
    /* We are done: */
    pswr_close_stream(fps);
  }

void draw_mst_edges(PSStream *fps, int N, r2_t st[], int P[], double C[])
  {
    int i;
    for (i = 0; i < N; i++)
      { int j = P[i];
        assert((j >= 0) && (j < N));
        if (i != j)
          { r2_t *si = &(st[i]);
            r2_t *sj = &(st[j]);
            pswr_segment(fps, si->c[0], si->c[1], sj->c[0], sj->c[1]);
          }
        else
          { assert(C[i] == +INF); }
      }
  }

void draw_sites(PSStream *fps, int N, r2_t st[], int P[], double C[])
  {
    double rvile, rroot;
    bool_t labels;
    if (N <= 20) 
      { rvile = 1.75; rroot = 2.00; labels = TRUE; }
    else
      { rvile = 0.50; rroot = 1.00; labels = FALSE; }
    if (labels)
      { pswr_set_label_font(fps, "CourierBold", 6.0); }
    int i;
    for (i = 0; i < N; i++) 
      { double x = st[i].c[0];
        double y = st[i].c[1];
        double r;
        if (C[i] == +INF)
          { pswr_set_fill_color(fps, 1.00,0.50,1.00); r = rroot; }
        else
          { pswr_set_fill_color(fps, 1.00,1.00,0.50); r = rvile; }
        pswr_dot(fps, x, y, r,  TRUE, TRUE); 
        if (labels)
          { char *lab = NULL;
            asprintf(&lab, "%d", i);
            pswr_label(fps, lab, x, y, 0.0, 0.5, 0.5);
            free(lab);
          }
      }
  }
  
double max_site_radius(int N, r2_t st[])
  { double r2max = 0.0;
    int i;
    for (i = 0; i < N; i++) 
      { double x = st[i].c[0];
        double y = st[i].c[1];
        double r2i = x*x + y*y;
        if (r2i > r2max) { r2max = r2i; }
      }
    return sqrt(r2max);
  }
  
PSStream *open_ps_stream(double r, char *prefix, bool_t eps_format)
  { 
    double mm = (72.0/25.4); /* One mm in pt. */
    double xfigsz = 150.00*mm; /* Figure X size excluding margin (pt). */
    double yfigsz = 150.00*mm; /* Figure Y size excluding margin (pt). */
    double fmrg = 3.0; /* Figure margin width (pt). */
    double pmrg = 2.0; /* Picture margin width (pt). */

    /* Add caption only if there is a user caption, or it is not EPS_FORMAT. */
    /* Select a good figure size: */
    PSStream *fps = pswr_new_stream
      ( /* prefix */                txtcat(prefix, "-"),
        /* file */                  NULL,
        /* eps_format */                   eps_format,
        /* docName */               "doc",
        /* paperSize */             "letter",
        /* landscape */             FALSE,
        /* hPageSize, vPageSize */  xfigsz + 2*fmrg, yfigsz + 2*fmrg
      );
    pswr_set_canvas_layout
      ( fps,
        /* hPicSize, vPicSize */     xfigsz, yfigsz,
        /* adjustPicSize */          FALSE,
        /* hPicMargin,vPicMargin */  pmrg, pmrg,
        /* captionLines */           (eps_format ? 0 : 1),  
        /* vCount, hCount */         0, 0  /* Let {pswr} choose it. */
      ); 
    return fps;
  }  
