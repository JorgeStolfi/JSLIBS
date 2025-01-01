/* see cpk_io.h */
/* Last edited on 2025-01-01 03:21:34 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <interval.h>
#include <r2.h>
#include <frgb.h>
#include <epswr.h>
#include <fget.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>

#include <cpk_io.h>
#include <cpk_main.h>
#include <cpk_coords.h>
#include <cpk_basic.h>
#include <hxg_eps.h>
#include <cpk_eps.h>

void cpk_plot_solution
  ( char *outDir,
    char *solTag,
    char *demTag,
    interval_t B[],
    cpk_domain_t *C,
    cpk_policy_t *P,
    r2_vec_t *V,
    double_vec_t *W,
    uint32_vec_t *J
  )
  {
    /* Compose output file name: */
    char *fileName = jsprintf("%s/%s-%s", outDir, solTag, demTag);
    
    /* Plot the domain and the candidates: */
    fprintf(stderr, "Writing Postscript file \"%s.eps\" ...\n", fileName);
    epswr_figure_t *eps = hxg_eps_new_figure(B, fileName);
    
    /* Show the domain: */
    cpk_plot_domain(eps, C, P);
    
    /* Show the demand points and their adequacy radii: */
    cpk_plot_demand(eps, &(C->Dem), P->rAuc);

    /* Show candidate points: */
    cpk_plot_candidates(eps, V, W);
    
    /* Show proposed auctions: */
    cpk_plot_proposed_auctions(eps, J, V, P->rAuc, P->dMin/2);
    
    /* Show statistics and such: */
    cpk_plot_solution_attributes(eps, B, P, W, J);

    /* A frame: */
    epswr_set_pen(eps, 0.8,0.5,0.2, 0.10, 0,0);
    epswr_frame(eps);
    
    epswr_end_figure(eps);
  }
  
r2_vec_t cpk_read_points(char *inDir, char *tag)
  {
    FILE *rd;
    { char *path = jsprintf("%s/%s.txt", inDir, tag);
      rd = open_read(path, TRUE);
      free(path);
    }
    
    /* Read file. */
    r2_vec_t P = r2_vec_new(100);
    uint32_t nP = 0;
    while(1)
      { /* Skip blanks on same line:: */
        fget_skip_spaces(rd);
        /* Check next character: */
        int32_t c = fgetc(rd);
        if (c == EOF)
          { break; }
        else if (c == '\n')
          { /* Blank line, ignore: */ }
        else if (c == '\015')
          { /* CR: ignore, like LF, but eat next LF if any: */
            c = fgetc(rd);
            if (c == EOF) { break; }
            if (c != '\n') { ungetc(c, rd); }
          }
        else if 
          ( (c == '+') || (c == '-') || (c == '.') || 
            ((c >= '0') && (c <= '9'))
          )
          { /* Number, must be a point: */
            ungetc(c, rd);
            double Lon = fget_double(rd);
            double Lat = fget_double(rd);
            r2_t p = (r2_t){{Lon,Lat}};
            /* Rest of line must be empty: */
            fget_eol(rd);
            /* Save point in list: */
            r2_vec_expand(&P, (vec_index_t)nP);
            P.e[nP] = p; 
            nP++;
          }
        else 
          { /* Non-number, non-blank: assume it is a comment, ignore: */
            do { c = fgetc(rd); } while ((c != '\n') && (c != EOF));
            if (c == EOF) { break; }
          }
      }
    /* Trim: */
    r2_vec_trim(&P, nP);
    return P;
  }

r2_vec_t cpk_read_polys(char *inDir, char *tag, bool_t closed)
  {
    FILE *rd;
    { char *path = jsprintf("%s/%s.txt", inDir, tag);
      rd = open_read(path, TRUE);
      free(path);
    }
    
    /* Read file. */
    r2_vec_t P = r2_vec_new(100);
    uint32_t nP = 0;
    int32_t fP = -1; /* First point of current component, or -1 outside components. */
    while(1)
      { /* Skip blanks on same line:: */
        fget_skip_spaces(rd);
        /* Check next character: */
        int32_t c = fgetc(rd);
        if (c == EOF)
          { break; }
        else if (c == '\n')
          { /* Blank line, ignore: */ }
        else if (c == '\015')
          { /* CR: ignore, like LF, but eat next LF if any: */
            c = fgetc(rd);
            if (c == EOF) { break; }
            if (c != '\n') { ungetc(c, rd); }
          }
        else if 
          ( (c == '+') || (c == '-') || (c == '.') || 
            ((c >= '0') && (c <= '9'))
          )
          { /* Number, must be a point: */
            ungetc(c, rd);
            double Lon = fget_double(rd);
            double Lat = fget_double(rd);
            r2_t p = (r2_t){{Lon,Lat}};
            /* Rest of line must be empty: */
            fget_eol(rd);
            /* Check component status: */
            if (fP == -1)
              { /* First point of a component. */
                if (nP > 0)
                  { /* Not the first component, must insert {{#}-vertex: */
                    r2_vec_expand(&P, (vec_index_t)nP);
                    P.e[nP] = (r2_t){{INF,INF}}; 
                    nP++;
                  }
                /* Remember this point: */
                fP = (int32_t)nP;
              }
            else if (closed)
              { /* Check for component closure: */
                if (r2_eq(&p, &(P.e[fP])))
                  { /* End of component: */
                    fP = -1; 
                  }
              }
            /* Save point in list: */
            r2_vec_expand(&P, (vec_index_t)nP);
            P.e[nP] = p; 
            nP++;
          }
        else 
          { /* Non-number, non-blank: assume it is a component separator: */
            if ((fP != -1) && closed) 
              { demand(FALSE, "component is not closed"); }
            fP = -1;
            do { c = fgetc(rd); } while ((c != '\n') && (c != EOF));
            if (c == EOF) { break; }
          }
      }
    /* End of input, check for proper closure: */
    if ((fP != -1) && closed)
      { demand(FALSE, "component is not closed"); }
    /* Trim: */
    r2_vec_trim(&P, nP);
    return P;
  }
  
void cpk_pick_ref_coords(r2_vec_t *P, double *refLon, double *refY)
  { 
    /* Simple-minded: take the averages over all finite vertices.
      Should be OK for any region in Brazil. For international 
      use we should check for wrap-around in longitude
      and/or in Northing. Some day ...
    */
    /* Compute the average longitude: */
    double sLon = 0;
    uint32_t nP = 0;
    for (int32_t i = 0; i < P->ne; i++) 
      { r2_t *pi = &(P->e[i]); 
        if (r2_is_finite(pi))
          { sLon += X(*pi); nP++; } 
      }
    /* If there are no finite points, return (0,0): */
    if (nP <= 0) { (*refLon) = 0; (*refY) = 0; return; }
    /* The reference longitude is the mean vertex longitude. */
    (*refLon) = sLon/nP;
    /* Map all points to EUTM, and compute the average Y: */
    double sY = 0; 
    for (int32_t i = 0; i < P->ne; i++) 
      { r2_t *Pi = &(P->e[i]); 
        if (r2_is_finite(Pi))
          { double Lon = X(*Pi);
            double Lat = Y(*Pi);
            double x, y;
            LL_to_EUTM(Lat, Lon, &x, &y, *refLon);
            sY += y;
          } 
      }
    (*refY) = sY/nP;
  }

r2_vec_t cpk_LL_to_EUTM(r2_vec_t *P, double refLon, double refY, double magnify)
  { r2_vec_t Q = r2_vec_new(P->ne);
    for (int32_t i = 0; i < P->ne; i++) 
      { r2_t *Pi = &(P->e[i]); 
        r2_t *Qi = &(Q.e[i]); 
        if (r2_is_finite(Pi))
          { double Lon = X(*Pi);
            double Lat = Y(*Pi);
            double x, y;
            LL_to_EUTM(Lat, Lon, &x, &y, refLon);
            (*Qi) = (r2_t){{x*magnify, (y - refY)*magnify}};
          }  
        else
          { (*Qi) = (*Pi); }
      }
    return Q;
  }

r2_vec_t cpk_EUTM_to_LL(r2_vec_t *P, double refLon, double refY, double magnify)
  { r2_vec_t Q = r2_vec_new(P->ne);
    for (int32_t i = 0; i < P->ne; i++) 
      { r2_t *Pi = &(P->e[i]); 
        r2_t *Qi = &(Q.e[i]); 
        if (r2_is_finite(Pi))
          { double x = X(*Pi)/magnify;
            double y = Y(*Pi)/magnify;
            double Lat, Lon;
            EUTM_to_LL(x, y + refY, &Lat, &Lon, refLon);
            (*Qi) = (r2_t){{Lon,Lat}};
          } 
        else
          { (*Qi) = (*Pi); }
      }
    return Q;
  }

void cpk_plot_domain(epswr_figure_t *eps, cpk_domain_t *C, cpk_policy_t *P)
  {
    /* Paint the polygon: */
    epswr_set_fill_color(eps, 1.00,0.95,0.93);
    cpk_eps_fill_polygon(eps, &(C->Urb));

    /* Unpaint the forbidden zones around existing/planned stations: */
    epswr_set_fill_color(eps, 1.00,1.00,1.00);
    cpk_eps_fill_circles(eps, &(C->Exs), NULL, P->dMin + P->rAuc);
    
    /* Unpaint the forbidden zones around predetermined auction centers: */
    epswr_set_fill_color(eps, 1.00,1.00,1.00);
    cpk_eps_fill_circles(eps, &(C->Auc), NULL, P->rAuc + P->dMin + P->rAuc);
    
    /* Unpaint the forbidden zone near intermunicipal borders: */
    cpk_eps_fill_sausage(eps, &(C->Mun), P->dMun);

    /* Unpaint the forbidden zone near international borders: */
    cpk_eps_fill_sausage(eps, &(C->Nat), P->dNat);

    /* Draw the polygon outline: */
    epswr_set_pen(eps, 0.50,0.25,0.00, 0.10, 0,0);
    cpk_eps_draw_polyline(eps, &(C->Urb));
    
    /* Draw the intermunicipal borders: */
    epswr_set_pen(eps, 0.50,0.00,1.00, 0.20, 0,0);
    cpk_eps_draw_polyline(eps, &(C->Mun));

    /* Draw the international borders: */
    epswr_set_pen(eps, 1.00,0.00,0.50, 0.30, 0,0);
    cpk_eps_draw_polyline(eps, &(C->Nat));

    /* Draw the exs./pln. stations: */
    frgb_t colorS = (frgb_t){{ 0.00f, 0.30f, 1.00f }};
    cpk_plot_stations(eps, &(C->Exs), 0.0, P->dMin/2, &colorS);

    /* Draw the predetermined auction points: */
    frgb_t colorA = (frgb_t){{ 0.50f, 0.30f, 0.00f }};
    cpk_plot_stations(eps, &(C->Auc), P->rAuc, P->dMin/2, &colorA);
  }

void cpk_plot_stations
  ( epswr_figure_t *eps, 
    r2_vec_t *P, 
    double rUnc, 
    double rRad,
    frgb_t *color
  )
  { 
    frgb_t z = (*color);
    
    /* Uncertainty circles (solid): */
    epswr_set_pen(eps, z.c[0],z.c[1],z.c[2], 0.10, 0,0);
    cpk_eps_draw_circles(eps, P, NULL, rUnc); 
    
    /* Nominal broadcast coverage areas (dashed): */
    epswr_set_pen(eps, z.c[0],z.c[1],z.c[2], 0.10, 1.0,0.5);
    cpk_eps_draw_circles(eps, P, NULL, rUnc+rRad); /* . */
    
    /* Station positions or auction centers (dots): */
    epswr_set_fill_color(eps, z.c[0],z.c[1],z.c[2]);
    cpk_eps_fill_dots(eps, P, 0.3);
    
    if (P->ne < 1000)
      { /* Write station indices: */
        r2_t disp = (r2_t){{ rRad/8, 0 }};
        epswr_set_pen(eps, z.c[0],z.c[1],z.c[2], 0.10, 0,0);
        cpk_eps_show_labels(eps, P, disp, 0.0,0.5, 8);
      }
  }

void cpk_plot_demand(epswr_figure_t *eps, r2_vec_t *P, double rDem)
  { 
    frgb_t z = (frgb_t){{ 0.00f, 0.60f, 0.00f }};
    
    /* Circle of adequate auction centers: */
    /* epswr_set_pen(eps, z.c[0],z.c[1],z.c[2], 0.10, 0,0); */
    /* cpk_eps_draw_circles(eps, P, NULL, rDem);  */
    
    /* Points of demand: */
    epswr_set_fill_color(eps, z.c[0],z.c[1],z.c[2]);
    cpk_eps_fill_dots(eps, P, 0.5);
  }
  
void cpk_plot_candidates(epswr_figure_t *eps, r2_vec_t *V, double_vec_t *W)
  {
    frgb_t z = (frgb_t){{ 0.70f, 0.40f, 0.00f }};
    
    uint32_t nV = V->ne;
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    for (int32_t i = 0; i < nV; i++)
      { r2_t *Vi = &(V->e[i]);
        double Ni = (i < W->ne ? W->e[i] : 0);
        double r = 0.10 + 0.05*floor(log(Ni+1)/log(2));
        if (Ni == 0)
          { epswr_set_fill_color(eps, 0.0,0.0,0.0); }
        else
          { epswr_set_fill_color(eps, z.c[0],z.c[1],z.c[2]); }
        epswr_dot(eps, X(*Vi), Y(*Vi), r, TRUE, FALSE);
      }
  }

void cpk_plot_proposed_auctions
  ( epswr_figure_t *eps, 
    uint32_vec_t *J, 
    r2_vec_t *V, 
    double rAuc, 
    double rRad
  )
  {
    uint32_t nV = V->ne;
    uint32_t nJ = J->ne;
    
    /* Extract the auction centers {U} with uncertainties {hU}: */
    uint32_t nU = nJ;
    r2_vec_t U = r2_vec_new(nU);
    { for (int32_t j = 0; j < nJ; j++)
        { uint32_t i = J->e[j]; 
          assert(i >= 0);
          assert(i < nV);
          U.e[j] = V->e[i];
        }
    }
    
    /* Now plot them: */
    frgb_t colorU = (frgb_t){{ 1.00f, 0.30f, 0.00f }};
    cpk_plot_stations(eps, &U, rAuc, rRad, &colorU);
    free(U.e);
  }

void cpk_plot_solution_attributes
  ( epswr_figure_t *eps, 
    interval_t B[], 
    cpk_policy_t *P,
    double_vec_t *W,
    uint32_vec_t *J
  )
  {
    uint32_t nJ = J->ne;
    
    /* Set font and choose position for key: */
    double hMin = 2.0*epswr_pt_per_mm;
    double hMax = hMin + 100*epswr_pt_per_mm;
    double vMin = 2.0*epswr_pt_per_mm;
    double vMax = vMin + 50*epswr_pt_per_mm;
    
    epswr_set_text_geometry(eps, FALSE, hMin, hMax, vMin, vMax, 0.0);
    epswr_set_pen(eps, 0.0,0.0,0.5, 0.10, 0.0, 0.0);
    double ptSize = 10;
    epswr_set_text_font(eps, "Courier", ptSize);
    
    /* Print the solution size: */
    { char *msg = jsprintf("Num auctions = %d", nJ);
      epswr_text(eps, msg, FALSE, 0.0, TRUE, FALSE);
      free(msg);
    }
    
    /* Print the solution weight: */
    { double WJ = 0.0;
      for (int32_t i = 0; i < nJ; i++) { WJ += W->e[J->e[i]]; }
      char *msg = jsprintf("Total weight = " WT_FFMT, WJ);
      epswr_text(eps, msg, FALSE, 0.0, TRUE, FALSE);
      free(msg);
    }
  }

void cpk_ui2_print(FILE *wr, ui2_t *x, char *fmt)
  { fputs("(", wr);
    for (int32_t i = 0; i < 2; i++)
      { if (i > 0) { fputs(" ", wr); }
        fprintf(wr, fmt, x->c[i]);
      }
    fputs(")", wr);
  }

