/* See pnmcolorfield.h
**
** Last edited on 2017-06-20 20:46:37 by stolfilocal
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
**
*/

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <jspnm.h>
#include <uint16_image.h>

#include <affirm.h>
#include <FBas.h>
#include <rn.h>

#include <colorfield.h>

/* INTERNAL PROTOTYPES */

ColorField *extract_coeffs
  ( int nb, 
    int chns, 
    double *W, 
    double tol
  );
  /* Assumes that {W[0..nb-1,0..chns-1]} contains the coefficients of
    an approximant relative to some function basis {bas[0..nb-1]}.
    Returns a {ColorField} {cf} comprising those terms whose
    coefficients, as determined by {A}, are significant, i.e. have
    some component greater than or equal to {tol} in absolute value.
    The {type}, {cols}, {rows}, and the {additive} flag of {cf} are
    left undefined. */
    
ColorField *make_uniform_color_field
  ( FBas *bas, 
    int chns,
    int cols, 
    int rows, 
    Triplet color,
    double tol
  );
  /* Returns a {ColorField} {cf} that evaluates everywhere to the
    given {color}. If the {color} is not significant, i.e. all its
    componnets are less than {tol} in absolute value, uses {Black}
    istead. The {type}, {cols}, {rows}, and the {additive} flag of
    {cf} are left undefined. */

int term_is_significant(Triplet *color, int chns, double tol);
  /* TRUE iff any component {color[0..chns-1]} is greater 
    than or equal to {tol} in absolute value. */

/* IMPLEMENTATIONS */

int cfld_debug = TRUE;

ColorField *compute_black_field
  ( ColorData *cd,    /* User specs of color field. */
    FBas *bas,        /* The approximation basis. */ 
    int chns,         /* Color channels to use. */
    int cols,         /* Number of columns in image. */
    int rows,         /* Number of rows in image. */
    double tol,       /* Discard terms that are less than this. */
    Triplet *inGamma  /* Gamma of color samples, per channel. */
  )
  { int np = cd->np;
    PixPos *pt = cd->pt;
    Triplet *raw = cd->color;
    
    /* Compute gamma-corrected black samples {cor[0..np-1]}: */
    Triplet cor[np];
    int k;
    for (k = 0; k < np; k++)
      { int ch;
        for (ch = 0; ch < chns; ch++)
          { cor[k].c[ch] = undo_in_gamma(raw[k].c[ch], inGamma->c[ch]); }
      }
      
    /* Fit color field: */
    ColorField *cf = fit_color_field(bas, np, pt, cor, chns, cols, rows, tol);
    cf->additive = TRUE;
    return cf;
  }
  
ColorField *compute_white_field
  ( ColorData *cd,     /* User specs of color field. */
    FBas *bas,         /* The approximation basis. */ 
    int chns,          /* Color channels to use. */
    int cols,          /* Number of columns in image. */
    int rows,          /* Number of rows in image. */
    double tol,        /* Discard terms that are less than this. */
    Triplet *inGamma,  /* Gamma of color samples, per channel. */
    ColorField *black, /* A previously computed {black} color field. */
    double minw        /* Minimum intensity for {white-black}. */
  )
  { int np = cd->np;
    PixPos *pt = cd->pt;
    Triplet *raw = cd->color;
    
    /* Compute gamma-corrected white samples {cor[0..np-1]}: */
    Triplet cor[np];
    int k;
    for (k = 0; k < np; k++)
      { int ch;
        int col = pt[k].c[0];
        int row = pt[k].c[1];
        Triplet bk = eval_color_field(black, col, row, chns);
        for (ch = 0; ch < chns; ch++)
          { double bkc = bk.c[ch];
            double whc = undo_in_gamma(raw[k].c[ch], inGamma->c[ch]);
            if (whc <= bkc) 
              { pm_error
                  ("white leq black [%d,%d] = (%.3f %.3f)", col,row,bkc,whc);
              }
            whc = whc - bkc;
            cor[k].c[ch] = log(whc < minw ? minw : whc);
          }
      }
      
    /* Fit color field: */
    ColorField *cf = fit_color_field(bas, np, pt, cor, chns, cols, rows, tol);
    cf->additive = FALSE;
    return cf;
  }  
    
ColorField *fit_color_field
  ( FBas *bas,      /* Basis to use. */  
    int np,         /* Number of samples. */
    PixPos *pt,     /* Sample positions. */ 
    Triplet *color, /* Sample values. */  
    int chns,       /* Color channels to use. */
    int cols,       /* Columns in image. */
    int rows,       /* Rows in image. */
    double tol      /* Discard terms that are less than this. */
  )
  { ColorField *cf;
    if (np == 1) 
      { /* Trivial (uniform) color field, no sweat: */
        cf = make_uniform_color_field(bas, chns, cols, rows, color[0], tol);
      }
    else 
      { /* Non-uniform color field, must fit coeffs to data. */
        /* Reformat data as expected by {FBas_Fit}: */
        int npclip = (np <= 15 ? np : 15);
        int nb = bas->m->roundsz(bas, npclip, TRUE);
        double P[2*np], V[chns*np], W[chns*nb];
        int k;
        for (k = 0; k < np; k++)
          { P[2*k+0] = 2*((double)pt[k].c[0])/((double)cols) - 1;
            P[2*k+1] = 2*((double)pt[k].c[1])/((double)rows) - 1;
            int ch;
            for (ch = 0; ch < chns; ch++)
              { V[chns*k+ch] = color[k].c[ch]; }
          }
        /* Fit the data: */
        FBas_Fit(np, 2, P, chns, V, bas, nb, W);

        /* Collect relevant terms: */
        cf = extract_coeffs(nb, chns, W, tol);
      }
    cf->bas = bas;
    cf->cols = cols;
    cf->rows = rows;
    cf->additive = TRUE;
    if(cfld_debug)
      { fprintf(stderr, "-- color field coeffs --\n");
        int t;
        for (t = 0; t < cf->nt; t++)
          { fprintf(stderr, "wt[%02d] = ", cf->index[t]); 
            print_triplet(stderr, "( ", chns, &(cf->wt[t]), " )");
            fprintf(stderr, "\n");
          }
        fprintf(stderr, "-- testing field at sample points --\n");
        int k;
        for (k = 0; k < np; k++)
          { PixPos *ptk = &(pt[k]);
            int col = ptk->c[0], row = ptk->c[1];
            Triplet *crk = &(color[k]);
            Triplet evk = eval_color_field(cf, col, row, chns);
            fprintf(stderr, "pt[%02d] = ", k);  
            fprintf(stderr, "(%04d,%04d)", col, row); 
            rn_gen_print(stderr, chns, crk->c, "%5.3f", "  dt = ( ", " ", " )"); 
            rn_gen_print(stderr, chns, evk.c, "%5.3f", "  ev =( ", " ", " )"); 
            fprintf(stderr, "\n");
          }
      }
    return cf;
  }

ColorField *extract_coeffs
  ( int nb, 
    int chns, 
    double *W, 
    double tol
  )
  { /* Temporary storage for relevant basis elements: */
    int nt = 0;
    int index[nb]; /* Relevant element indices. */
    Triplet wt[nb]; /* Corresponding coefficients. */
    /* Gather significant terms: */
    int i;
    for (i = 0; i < nb; i++)
      { /* Repackage coeffs of {bas[i]} as a color triplet: */
        Triplet color = Black;
        double *Wi = &(W[chns*i]); 
        int ch;
        for (ch = 0; ch < chns; ch++) { color.c[ch] = Wi[ch]; }
        /* Check if significant: */
        if (term_is_significant(&color, chns, tol))
          { /* Store the coefficient: */
            index[nt] = i;  wt[nt] = color; nt++;
          }
      }
    /* Build colorfield structure: */
    ColorField *cf = (ColorField *)malloc(sizeof(ColorField));
    affirm(cf != NULL, "out of mem");
    cf->nt = nt;
    cf->index = (int *)malloc(nt * sizeof(index));
    affirm((nt == 0) || (cf->index != NULL), "out of mem");
    cf->wt = (Triplet *)malloc(nt * sizeof(Triplet));
    affirm((nt == 0) || (cf->wt != NULL), "out of mem");
    int k;
    for (k = 0; k < nt; k++) 
      { cf->index[k] = index[k]; cf->wt[k] = wt[k]; }
    return cf;
  }
 
ColorField *make_uniform_color_field
  ( FBas *bas, 
    int chns,
    int cols, 
    int rows, 
    Triplet color,
    double tol
  )
  { int nt = (term_is_significant(&color, chns, tol) ? 1 : 0);
    /* Build colorfield structure: */
    ColorField *cf = (ColorField *)malloc(sizeof(ColorField));
    affirm(cf != NULL, "out of mem");
    cf->nt = nt;
    cf->index = (int *)malloc(nt * sizeof(int));
    affirm((nt == 0) || (cf->index != NULL), "out of mem");
    cf->wt = (Triplet *)malloc(nt * sizeof(Triplet));
    affirm((nt == 0) || (cf->wt != NULL), "out of mem");
    /* Assumes that {bas[0](p) == 1} for all {p}: */
    if (nt > 0) 
      { cf->index[0] = 0; 
        int ch;
        for (ch = 0; ch < 3; ch++)
          { cf->wt[0].c[ch] = (ch < chns ? color.c[ch] : 0.0); }
      }
    return cf;
  }
  
Triplet eval_color_field(ColorField *cf, int col, int row, int chns)
  { Triplet color = Black;
    int nt = cf->nt, cols = cf->cols, rows = cf->rows;
    FBas *bas = cf->bas;
    int k, ch;
    double p[2];
    p[0] = 2*((double)col)/((double)cols) - 1.0;
    p[1] = 2*((double)row)/((double)rows) - 1.0;
    for (k = 0; k < nt; k++)
      { /* Get the index and coeff of relevant term {[k]}: */
        int idk = cf->index[k];
        Triplet *wk = &(cf->wt[k]);
        /* Evaluate that basis element at {(col,row)}: */
        double bvk = bas->m->eval(bas, idk, 2, p);
        /* fprintf(stderr, "basis[%d](%d,%d) = %12.6f\n", idk, col, row, bvk); */
        for (ch = 0; ch < chns; ch++)
          { color.c[ch] += bvk*wk->c[ch]; }
      }
    if (! cf->additive) 
      { for (ch = 0; ch < chns; ch++)
          { color.c[ch] = exp(color.c[ch]); }
      }
    return color;
  }

int eq_triplet(Triplet a, Triplet b, int chns)
  { int ch;
    for (ch = 0; ch < chns; ch++)
      { if (a.c[ch] != b.c[ch]) return FALSE; }
    return TRUE;
  }

int field_is_constant(ColorField *cf, Triplet color, int chns)
  { /* If there are no relevant terms, the field is uniformly black: */
    if (cf->nt == 0) { return eq_triplet(color, Black, chns); }
    /* If there is more than one term, the field is not uniform: */
    if (cf->nt > 1) { return FALSE; }
    /* If the only term is not {basis[0]}, the field is not uniform: */
    if (cf->index[0] != 0) { return FALSE; }
    /* Field is uniform, check the color: */
    Triplet *w0 = &(cf->wt[0]);
    return eq_triplet(*w0, color, chns);
  }

int term_is_significant(Triplet *color, int chns, double tol)
  { int ch;
    for (ch = 0; ch < chns; ch++) 
      { double v = fabs(color->c[ch]); 
        if (v >= tol) { return TRUE; }
      }
    return FALSE;
  }

void cfld_debug_triplet(char *label, int col, int row, int chns, Triplet *fv, char *tail)
  { 
    if (cfld_debug) 
      { fprintf(stderr, "%s[%3d][%3d] = ", label, row, col);
        print_triplet(stderr, "( ", chns, fv, " )");
        fprintf(stderr, "%s", tail);
      }
  }
       
      
