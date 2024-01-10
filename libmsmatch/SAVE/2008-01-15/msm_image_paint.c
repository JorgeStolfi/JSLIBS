/* See {msm_image_paint.h} */
/* Last edited on 2008-01-12 08:24:52 by stolfi */

#define msm_image_paint_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_seq_desc.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>
#include <msm_image.h>
#include <msm_dyn.h>
#include <msm_image_paint.h>

#include <vec.h>
#include <float_image.h>
#include <jsmath.h>
#include <affirm.h>
#include <assert.h>

#include <stdint.h>
#include <math.h>

void msm_image_seq_seq_score_paint
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp, 
    msm_rung_func_t *score, 
    msm_image_t *img,
    int xMin, 
    int yMin
  )
  { int NC = img->fim->sz[0];
    demand(NC == 1, "image must be monochromatic");
    /* Get the sequence sizes {nx,ny}: */
    int nx = xp->nbas;
    int ny = yp->nbas;
    /* Clear all totals: */
    float_image_fill(img->fim, 0.0);
    /* Value to add to pixel for each sample pair: */
    float wt = 1.0/(img->scale*img->scale);
    /* Count coincidence points: */
    int iy, ix;
    for (iy = 0; iy < ny; iy++)
      { for (ix = 0; ix < nx; ix++)
          { msm_rung_t g = (msm_rung_t){{ ix, iy }};
            double sc = score(&g);
            /* fprintf(stderr, "    score = %10.6f\n", sc); */
            msm_image_add(img, 0, ix-xMin, iy-yMin, wt*sc);
          }
      }
  }

void msm_image_pairing_paint
  ( msm_pairing_t *p,
    int nx,
    int ny,
    float v[],
    bool_t add,
    msm_image_t *img,
    int xMin, 
    int yMin
  )
  { /* Get number of channels in image: */
    int NC = img->fim->sz[0];
    /* Get number of rungs in candidate's pairing: */
    int ng = msm_pairing_num_rungs(p);
    /* Enumerate sample pairs in candidate: */
    int ig, ic;
    for (ig = 0; ig < ng; ig++)
      { msm_rung_t g = msm_pairing_get_rung(p, ig);
        /* Get indices into sequences: */
        int32_t ix, iy;
        ix = g.c[0]; iy = g.c[1];
        if (add) 
          { /* Assumes that at most {scale} rungs will fall into the same pixel: */
            for (ic = 0; ic < NC; ic++)
              { msm_image_add(img, ic, ix-xMin, iy-yMin, v[ic]/img->scale); }
          }
        else
          { /* Set whole pixel if any rung falls into it: */
            for (ic = 0; ic < NC; ic++)
              { msm_image_set(img, ic, ix-xMin, iy-yMin, v[ic]); }
          }
      }
  }
 
void msm_image_cand_paint
  ( msm_cand_t *cd, 
    float v[], 
    bool_t add, 
    msm_image_t *img, 
    int xMin, 
    int yMin
  )
  { /* Get sequence sizes {nx,ny}: */
    int nx = cd->seq[0].nbas;
    int ny = cd->seq[1].nbas;
    /* Paint rungs: */
    msm_image_pairing_paint(cd->pr, nx, ny, v, add, img, xMin, yMin);
 }

void msm_image_cand_vec_paint
  ( msm_cand_vec_t *cdv,
    msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp,
    msm_image_t *img,
    int xMin, 
    int yMin
  )
  {
    int NC = img->fim->sz[0];
    demand(NC == 1, "image must be monochromatic");
    /* Fill the image with "no value": */
    float_image_fill(img->fim, (float)(-INF));
    /* Enumerate candidates for {xp,yp} and accumulate the pairs: */
    int ic; 
    for (ic = 0; ic < cdv->nel; ic++)
      { /* Get the candidate {cd} with index {ic}: */
        msm_cand_t *cd = &(cdv->el[ic]);
        /* Compute average score {avsc} per rung: */
        int ng = msm_pairing_num_rungs(cd->pr);
        double avsc = cd->score/((double)ng);
        /* Compute the value to paint into {img}: */
        float v = ((double)avsc)/((double)img->scale);
        /* Check whether seqences have same size: */
        demand(cd->seq[0].nbas == xp->nbas, "inconsistent xp->nbas");
        demand(cd->seq[1].nbas == yp->nbas, "inconsistent yp->nbas");
        /* Paint it: */
        msm_image_cand_paint(cd, &v, TRUE, img, xMin, yMin);
      }
  }
  
void msm_image_dyn_tableau_scores_paint
  ( int nx,                /* Length of first sequence. */
    int ny,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau. */
    msm_image_t *img,
    int xMin, 
    int yMin
  )
  { float_image_fill(img->fim, -INF);

    /* Paint tableau {tb} into float image {img}: */
    int dr;
    for (dr = 0; dr < tb->nr; dr++)
      { /* Compute coordinnate {r} of rung: */
        int r = tb->rMin + dr;
        /* Compute range {sMin..sMax} of coord {s} for this {r}: */
        int sMin, sMax;
        msm_dyn_tableau_get_s_range(tb, r, &sMin, &sMax);
        if (sMin <= sMax)
          { demand((sMin + r) % 2 == 0, "{sMin} has wrong parity");
            demand((sMax + r) % 2 == 0, "{sMax} has wrong parity");
            /* Compute number {ns} of distinct coords {s} for this {r}: */
            int ns = (sMax - sMin)/2 + 1;
            /* Enumerate {s} coords for this {r}: */
            int ds;
            for (ds = 0; ds < ns; ds++)
              { /* Compute {s} coord: */
                int s = sMin + 2*ds;
                /* Compute index of entry in {tb->ev} and get its address {ep}: */
                int k = dr*tb->ns + ds;
                msm_dyn_entry_t *ep = &(tb->ev.el[k]);
                /* Compute the rung {g = (ix,iy)} corresponding to this entry: */ 
                assert ((r + s) % 2 == 0); /* Parity constraint. */
                int ix = (r + s)/2;
                int iy = (r - s)/2;
                /* Paint score into image (use {max} when reducing): */
                double v = ep->score;
                msm_image_max(img, 0, ix-xMin, iy-yMin, v);
              }
          }
      }
  }

void msm_image_dyn_tableau_pairing_paint
  ( int nx,                /* Length of first sequence. */
    int ny,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau. */
    msm_rung_t gopt,       /* Optimal end-rung or {msm_rung_none} */ 
    float clr[],
    msm_image_t *img,
    int xMin, 
    int yMin
  )
  {
    int NC = img->fim->sz[0];
    msm_rung_t g = gopt;
    while (! msm_rung_is_none(&g))
      { /* Get the indices {ix,iy} from the rung: */
        int ix = g.c[0];
        int iy = g.c[1];
        /* Paint the rung into {cim} with color {clr}: */
        int c;
        for (c = 0; c < NC; c++) 
          { msm_image_set(img, c, ix-xMin, iy-yMin, clr[c]); }
        msm_dyn_entry_t *ep = msm_dyn_tableau_get_entry_address(tb, g);
        affirm(ep != NULL, "optmum path falls off the tableau");
        g = ep->prev;
      }
  }

void msm_image_seq_periods_paint
  ( int nx,                /* Length of first sequence. */
    int ny,                /* Length of second sequence. */
    float clr[],
    msm_image_t *img,
    int xMin, 
    int yMin
  )
  { /* Get the virtual image dimensions {NVX,NVY}: */
    int NVX = img->NVX;
    int NVY = img->NVY;
    /* Get the float image dimensions {NFC,NFX,NFY}: */
    int NC = img->fim->sz[0];
    int NFX = img->fim->sz[1];
    int NFY = img->fim->sz[2];
    /* Check compatibility with {nx,ny}: */
    demand(NVX % nx == 0, "bad image X size");
    demand(NVY % ny == 0, "bad image Y size");
    /* Compute the first virtual X and Y coordinates {ix,iy}. */
    int xIni = iceil(xMin, nx); assert((xIni >= xMin) && (xIni < xMin + ny));
    int yIni = iceil(yMin, ny); assert((yIni >= yMin) && (yIni < yMin + ny));
    /* Compute the number of plottable grid lines {mx,my}. */
    int mx = NVX/nx;
    int my = NVY/ny;
    int k; 
    /* Plot the vertical lines: */
    for (k = 0; k < mx; k++)
      { /* Compute the virtual X coordinate of line number {k}: */
        int x = xIni + k*nx;
        /* Compute the float image X coordinate of line number {k}: */
        int ixp = (x - xMin)/img->scale;
        /* Paint all its pixels: */
        int iyp;
        for (iyp = 0; iyp < NFY; iyp++)
          { int c;
            for (c = 0; c < NC; c++)
              { float_image_set_sample(img->fim, c, ixp, iyp, clr[c]); }
          }
      }
    /* Plot the horizontal lines: */
    for (k = 0; k < my; k++)
      { /* Compute the virtual Y coordinate of line number {k}: */
        int y = yIni + k*ny;
        /* Compute the float image Y coordinate of line number {k}: */
        int iyp = (y - yMin)/img->scale;
        /* Paint all its pixels: */
        int ixp;
        for (ixp = 0; ixp < NFX; ixp++)
          { int c;
            for (c = 0; c < NC; c++)
              { float_image_set_sample(img->fim, c, ixp, iyp, clr[c]); }
          }
      }
  }

