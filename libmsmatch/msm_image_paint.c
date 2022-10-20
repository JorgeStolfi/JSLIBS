/* See {msm_image_paint.h} */
/* Last edited on 2022-10-20 06:38:50 by stolfi */

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
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1, 
    msm_rung_score_proc_t *score, 
    msm_image_t *img,
    int32_t i0Min, 
    int32_t i1Min
  )
  { int32_t NC = (int32_t)img->fim->sz[0];
    demand(NC == 1, "image must be monochromatic");
    /* Get the sequence sizes {n0,n1}: */
    int32_t n0 = seq0->size;
    int32_t n1 = seq1->size;
    /* Clear all totals: */
    float_image_fill(img->fim, 0.0);
    /* Value to add to pixel for each sample pair: */
    double wt = 1.0/(img->scale*img->scale);
    /* Count coincidence points: */
    int32_t i1, i0;
    bool_t debug = TRUE;
    if(debug){fprintf(stderr,"\n");}
    for (i1 = 0; i1 < n1; i1++)
      { for (i0 = 0; i0 < n0; i0++)
          { msm_rung_t g = (msm_rung_t){{ i0, i1 }};
            double sc = score(seq0, seq1, &g);
            /* fprintf(stderr, "    score = %10.6f\n", sc); */
            msm_image_add(img, 0, i0-i0Min, i1-i1Min, (float)(wt*sc));
          }
      }
  }

void msm_image_pairing_paint
  ( msm_pairing_t *p,
    int32_t n0,
    int32_t n1,
    float v[],
    bool_t add,
    msm_image_t *img,
    int32_t i0Min, 
    int32_t i1Min
  )
  { /* Get number of channels in image: */
    int32_t NC = (int32_t)img->fim->sz[0];
    /* Get number of rungs in candidate's pairing: */
    int32_t ng = msm_pairing_num_rungs(p);
    /* Enumerate sample pairs in candidate: */
    int32_t ig, ic;
    for (ig = 0; ig < ng; ig++)
      { msm_rung_t g = msm_pairing_get_rung(p, ig);
        /* Get indices into sequences: */
        int32_t i0, i1;
        i0 = g.c[0]; i1 = g.c[1];
        if (add) 
          { /* Assumes that at most {scale} rungs will fall into the same pixel: */
            for (ic = 0; ic < NC; ic++)
              { float vs = (float)(v[ic])/(float)(img->scale);
                msm_image_add(img, ic, i0-i0Min, i1-i1Min, vs);
              }
          }
        else
          { /* Set whole pixel if any rung falls into it: */
            for (ic = 0; ic < NC; ic++)
              { msm_image_set(img, ic, i0-i0Min, i1-i1Min, v[ic]); }
          }
      }
  }
 
void msm_image_cand_paint
  ( msm_cand_t *cd, 
    float v[], 
    bool_t add, 
    msm_image_t *img, 
    int32_t i0Min, 
    int32_t i1Min
  )
  { /* Get sequence sizes {n0,n1}: */
    int32_t n0 = cd->seq[0].size;
    int32_t n1 = cd->seq[1].size;
    /* Paint rungs: */
    msm_image_pairing_paint(cd->pr, n0, n1, v, add, img, i0Min, i1Min);
 }

void msm_image_cand_vec_paint
  ( msm_cand_vec_t *cdv,
    msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    msm_image_t *img,
    int32_t i0Min, 
    int32_t i1Min
  )
  {
    int32_t NC = (int32_t)img->fim->sz[0];
    demand(NC == 1, "image must be monochromatic");
    /* Fill the image with "no value": */
    float_image_fill(img->fim, (float)(-INF));
    /* Enumerate candidates for {seq0,seq1} and accumulate the pairs: */
    int32_t ic; 
    for (ic = 0; ic < cdv->ne; ic++)
      { /* Get the candidate {cd} with index {ic}: */
        msm_cand_t *cd = &(cdv->e[ic]);
        /* Compute average score {avsc} per rung: */
        int32_t ng = msm_pairing_num_rungs(cd->pr);
        double avsc = cd->score/((double)ng);
        /* Compute the value to paint into {img}: */
        float v = (float)(((double)avsc)/((double)img->scale));
        /* Check whether seqences have same size: */
        demand(cd->seq[0].size == seq0->size, "inconsistent seq0->size");
        demand(cd->seq[1].size == seq1->size, "inconsistent seq1->size");
        /* Paint it: */
        msm_image_cand_paint(cd, &v, TRUE, img, i0Min, i1Min);
      }
  }
  
void msm_image_dyn_tableau_scores_paint
  ( int32_t n0,                /* Length of first sequence. */
    int32_t n1,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau. */
    msm_image_t *img,
    int32_t i0Min, 
    int32_t i1Min
  )
  { float_image_fill(img->fim, -INF);

    /* Paint tableau {tb} into float image {img}: */
    int32_t dr;
    for (dr = 0; dr < tb->nr; dr++)
      { /* Compute coordinnate {r} of rung: */
        int32_t r = tb->rMin + dr;
        /* Compute range {sMin..sMax} of coord {s} for this {r}: */
        int32_t sMin, sMax;
        msm_dyn_tableau_get_s_range(tb, r, &sMin, &sMax);
        if (sMin <= sMax)
          { demand((sMin + r) % 2 == 0, "{sMin} has wrong parity");
            demand((sMax + r) % 2 == 0, "{sMax} has wrong parity");
            /* Compute number {ns} of distinct coords {s} for this {r}: */
            int32_t ns = (sMax - sMin)/2 + 1;
            /* Enumerate {s} coords for this {r}: */
            int32_t ds;
            for (ds = 0; ds < ns; ds++)
              { /* Compute {s} coord: */
                int32_t s = sMin + 2*ds;
                /* Compute index of entry in {tb->ev} and get its address {ep}: */
                int32_t k = dr*tb->ns + ds;
                msm_dyn_entry_t *ep = &(tb->ev.e[k]);
                /* Compute the rung {g = (i0,i1)} corresponding to this entry: */ 
                assert ((r + s) % 2 == 0); /* Parity constraint. */
                int32_t i0 = (r + s)/2;
                int32_t i1 = (r - s)/2;
                /* Paint score into image (use {max} when reducing): */
                double v = ep->score;
                msm_image_max(img, 0, i0-i0Min, i1-i1Min, (float)v);
              }
          }
      }
  }

void msm_image_dyn_tableau_pairing_paint
  ( int32_t n0,                /* Length of first sequence. */
    int32_t n1,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau. */
    msm_rung_t gopt,       /* Optimal end-rung or {msm_rung_none} */ 
    float clr[],
    msm_image_t *img,
    int32_t i0Min, 
    int32_t i1Min
  )
  {
    int32_t NC = (int32_t)img->fim->sz[0];
    msm_rung_t g = gopt;
    while (! msm_rung_is_none(&g))
      { /* Get the indices {i0,i1} from the rung: */
        int32_t i0 = g.c[0];
        int32_t i1 = g.c[1];
        /* Paint the rung into {cim} with color {clr}: */
        int32_t c;
        for (c = 0; c < NC; c++) 
          { msm_image_set(img, c, i0-i0Min, i1-i1Min, clr[c]); }
        msm_dyn_entry_t *ep = msm_dyn_tableau_get_entry_address(tb, g);
        affirm(ep != NULL, "optmum path falls off the tableau");
        g = ep->prev;
      }
  }

void msm_image_grid_paint
  ( int32_t step0,
    int32_t step1,
    float clr[],
    msm_image_t *img,
    int32_t i0Min, 
    int32_t i1Min
  )
  { /* Get the virtual image dimensions {NVX,NVY}: */
    int32_t NVX = img->NV[0];
    int32_t NVY = img->NV[1];
    /* Get the actual image dimensions {NFC,NFX,NFY}: */
    int32_t NC =  (int32_t)img->fim->sz[0];
    int32_t NFX = (int32_t)img->fim->sz[1];
    int32_t NFY = (int32_t)img->fim->sz[2];
    /* Compute the first virtual X and Y coordinates {i0,i1}. */
    int32_t xIni = (int32_t)iceil(i0Min, step0); assert((xIni >= i0Min) && (xIni < i0Min + step1));
    int32_t yIni = (int32_t)iceil(i1Min, step1); assert((yIni >= i1Min) && (yIni < i1Min + step1));
    /* Compute the number of plottable grid lines {mx,my}. */
    int32_t mx = NVX/step0;
    int32_t my = NVY/step1;
    int32_t k; 
    /* Plot the vertical lines: */
    for (k = 0; k < mx; k++)
      { /* Compute the virtual X coordinate of line number {k}: */
        int32_t x = xIni + k*step0;
        /* Compute the float image X coordinate of line number {k}: */
        int32_t ixp = (x - i0Min)/img->scale;
        /* Paint all its pixels: */
        int32_t iyp;
        for (iyp = 0; iyp < NFY; iyp++)
          { int32_t c;
            for (c = 0; c < NC; c++)
              { float_image_set_sample(img->fim, c, ixp, iyp, clr[c]); }
          }
      }
    /* Plot the horizontal lines: */
    for (k = 0; k < my; k++)
      { /* Compute the virtual Y coordinate of line number {k}: */
        int32_t y = yIni + k*step1;
        /* Compute the float image Y coordinate of line number {k}: */
        int32_t iyp = (y - i1Min)/img->scale;
        /* Paint all its pixels: */
        int32_t ixp;
        for (ixp = 0; ixp < NFX; ixp++)
          { int32_t c;
            for (c = 0; c < NC; c++)
              { float_image_set_sample(img->fim, c, ixp, iyp, clr[c]); }
          }
      }
  }

