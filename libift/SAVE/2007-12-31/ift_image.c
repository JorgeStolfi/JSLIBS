/* Implementation of ift_image.h */
/* Last edited on 2007-12-26 20:52:46 by stolfi */

#include "ift_image.h"
#include <jspnm_image.h>
#include <math.h>

#include <jspnm.h>
#include <jspnm_image.h>

void ift_check_maxval(unsigned int maxval, unsigned int maxmaxval);
  /* Bombs out if {maxval} is zero or greater than {maxmaxval}. */

void ift_check_maxval(unsigned int maxval, unsigned int maxmaxval)
  { if ((maxval < 1) || (maxval > maxmaxval)) 
      { fprintf(stderr, "maxval = %d\n", maxval); 
        IFT_ERROR("bad maxval");
      }
  }

void ift_set_values_from_image(pnm_image_t *img, ImageGraph *G)
  { int col, row, c; 
    PixelNode *pg = &(G->node[0]);
    PixelValue yy;
    double fmaxval = (double)(img->maxval);
    if (img->cols != G->cols) { IFT_ERROR("cols mismatch"); }
    if (img->rows != G->rows) { IFT_ERROR("rows mismatch"); }
    G->channels = img->chns;
    /* Initialize unused channels: */
    for (c = G->channels; c < MAX_CHANNELS; c++) { yy.c[c] = 0; }
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { pnm_sample_t *pm = &(img->smp[row][col*img->chns]);
          for (c = 0; c < img->chns; c++) { yy.c[c] = ((double)(*pm))/fmaxval; pm++; }
          pg->y = yy;
          pg++;
        }
  }

void ift_set_seeds_from_image(pnm_image_t *seed_img, ImageGraph *G)
  { int col, row;
    PixelNode *pg = &(G->node[0]);
    if (seed_img->chns != 1) { IFT_ERROR("seed image must be monochromatic"); }
    if (seed_img->cols != G->cols) { IFT_ERROR("cols mismatch"); }
    if (seed_img->rows != G->rows) { IFT_ERROR("rows mismatch"); }
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { unsigned int pm = seed_img->smp[row][col];
          if (pm > MAX_LABEL) 
            { IFT_ERROR("seed label too big"); }
          else
            { pg->L = pm; }
          pg++;
        }
  }

pnm_image_t *ift_get_cost_image(ImageGraph *G, PathCost maxcost, pnm_sample_t maxval)
  { int col, row;
    double scale;
    PixelNode *pg = &(G->node[0]);
    pnm_image_t *img = pnm_image_new(G->cols, G->rows, 1);
    ift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);
    img->maxval = maxval;
    if (maxcost < 1.0)  { IFT_ERROR("bad maxcost"); }
    scale = ((double)(maxval - 1))/((double)maxcost);
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { pnm_sample_t *pm = &(img->smp[row][col]);
          if (pg->C == INFINITE_PATH_COST)
            { (*pm) = img->maxval; }
          else
            { double y = scale*pg->C;
              if (y > (double)(img->maxval-1)) 
                { (*pm) = img->maxval-1; }
              else
                { (*pm) = (int)(y + 0.5); }
            }
          pg++;
        }
    return img;
  }

pnm_image_t *ift_get_pred_image(ImageGraph *G)
  { int col, row, dcol, drow;
    PixelNode *pg = &(G->node[0]);
    pnm_image_t *img = pnm_image_new(G->cols, G->rows, 1);
    img->maxval = PNM_FILE_MAX_MAXVAL;
    if (img->maxval < 65535) { IFT_ERROR("maxval too small"); }
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { pnm_sample_t *pm = &(img->smp[row][col]);
          PixelNode *pp = (PixelNode *)pg->P;
          if (pp == NULL)
            { dcol = 0; drow = 0; }
          else
            { dcol = pp->col - col; drow = pp->row - row; }
          if ((dcol < -127) || (dcol > 127)) { IFT_ERROR("pred dcol out of range"); }
          if ((drow < -127) || (drow > 127)) { IFT_ERROR("pred drow out of range"); }
          (*pm) = ((drow + 128) << 8) | (dcol + 128);
          pg++;
        }
    return img;
  }

pnm_image_t *ift_get_label_image(ImageGraph *G, pnm_sample_t maxval)
  { int col, row;
    PixelNode *pg = &(G->node[0]);
    pnm_image_t *img = pnm_image_new(G->cols, G->rows, 1);
    ift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);
    img->maxval = maxval;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { pnm_sample_t *pm = &(img->smp[row][col]);
          if (pg->L > maxval) { IFT_ERROR("label too big"); }
          if (pg->R == NULL) { IFT_ERROR("null root"); }
          if (pg->L != (pg->R)->L) { IFT_ERROR("inconsistent label"); }
          (*pm) = pg->L;
          pg++;
        }
    return img;
  }

pnm_image_t *ift_get_root_image(ImageGraph *G)
  { int col, row;
    PixelNode *pg = &(G->node[0]);
    pnm_image_t *img = pnm_image_new(G->cols, G->rows, 1);
    img->maxval = 2;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { pnm_sample_t *pm = &(img->smp[row][col]);
          if ((pg->P == NULL) != (pg->R == pg)) { IFT_ERROR("P/R inconsistency"); }
          (*pm) = (pg->P == NULL ? 2 : 0);
          pg++;
        }
    return img;
  }

pnm_image_t *ift_get_spread_image(ImageGraph *G, pnm_sample_t maxval)
  { pnm_image_t *img;
    int col, row, ch;
    PixelNode *pg = &(G->node[0]);
    double fmaxval = (double)maxval;
    /* Choose imag etype: */
    ift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);

    /* Allocate image and fill it out: */
    img = pnm_image_new(G->cols, G->rows, G->channels);
    img->maxval = maxval;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { PixelNode *pr = pg->R;
          if (pr == NULL) { IFT_ERROR("R is null"); }
          if ((pr == pg) != (pg->P == NULL)) { IFT_ERROR("P/R inconsistency"); }
          pnm_sample_t *pm = &(img->smp[row][col*img->chns]);
          for (ch = 0; ch < G->channels; ch++)
            { double yy = pr->y.c[ch];
              int pv;
              if (yy <= 0.0) 
                { pv = 0; }
              else if (yy >= 1.0)
                { pv = maxval; }
              else
                { pv = (int)floor(0.5 + yy*fmaxval); }
              (*pm) = pv;
	      pm++;
            }
          pg++;
        }
    return img;
  }

pnm_image_t *ift_get_single_label_image
  ( ImageGraph *G, 
    SeedLabel label, 
    pnm_sample_t bg[], 
    pnm_sample_t maxval
  )
  { pnm_image_t *img;
    int col, row, ch;
    PixelNode *pg = &(G->node[0]);
    double fmaxval = (double)maxval;
    /* Choose imag etype: */
    ift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);

    /* Allocate image and fill it out: */
    img = pnm_image_new(G->cols, G->rows, G->channels);
    img->maxval = maxval;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { SeedLabel pglab = pg->L;
          pnm_sample_t *pm = &(img->smp[row][col*img->chns]);
          for (ch = 0; ch < G->channels; ch++)
            { 
              int pv;
              if (pglab == label) 
                { double yy = pg->y.c[ch];
                  if (yy <= 0.0) 
                    { pv = 0; }
                  else if (yy >= 1.0)
                    { pv = maxval; }
                  else
                    { pv = (int)floor(0.5 + yy*fmaxval); }
                }
              else
                { pv = bg[ch]; }
              (*pm) = pv;
	      pm++;
            }
          pg++;
        }
    return img;
  }

void ift_write_boxes(FILE *wr, ImageGraph *G, pnm_sample_t maxval, int margin)
  { int col, row;
    
    /* Max and min pixel indices for each label: */
    int *colmin = (int *)malloc((maxval+1)*sizeof(int));
    int *colmax = (int *)malloc((maxval+1)*sizeof(int));
    int *rowmin = (int *)malloc((maxval+1)*sizeof(int));
    int *rowmax = (int *)malloc((maxval+1)*sizeof(int));
    PixelNode *pg = &(G->node[0]);
    
    ift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);

    /* Initialize max/min accumulators: */
    if ((colmin == NULL) || (colmax == NULL)) { IFT_ERROR("out of mem"); }
    if ((rowmin == NULL) || (rowmax == NULL)) { IFT_ERROR("out of mem"); }
    { int lab;
      for (lab = 0; lab <= maxval; lab++) 
        { colmin[lab] = rowmin[lab] = 1000000; 
          colmax[lab] = rowmax[lab] = 0;
        }
    }
      
    /* Scan image and collect max/min indices per label: */
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { SeedLabel lab = pg->L;
          if (pg->L > maxval) { IFT_ERROR("label too big"); }
          if (pg->R == NULL) { IFT_ERROR("null root"); }
          if (pg->L != (pg->R)->L) { IFT_ERROR("inconsistent label"); }
          if (col < colmin[lab]) { colmin[lab] = col; }
          if (col > colmax[lab]) { colmax[lab] = col; }
          if (row < rowmin[lab]) { rowmin[lab] = row; }
          if (row > rowmax[lab]) { rowmax[lab] = row; }
          pg++;
        }
      
    /* Print region boxes: */
    { int lab;
      for (lab = 0; lab <= maxval; lab++) 
        { int colLO = colmin[lab], colHI = colmax[lab] + 1;
          int rowLO = rowmin[lab], rowHI = rowmax[lab] + 1;
          fprintf(wr, "%5d ", lab);
          /* Apply margin expansion/contraction: */
          if ((colLO < colHI) && (rowLO < rowHI))
            { colLO -= margin; colHI += margin;
              rowLO -= margin; rowHI += margin;
              if (colLO < 0) { colLO = 0; }
              if (colHI > G->cols) { colHI = G->cols; }
              if (rowLO < 0) { rowLO = 0; }
              if (rowHI > G->rows) { rowHI = G->rows; }
            }
          /* Normalize empty boxes: */
          if ((colLO >= colHI) || (rowLO >= rowHI))
            { colLO = colHI = 0; rowLO = rowHI = 0; }
          /* Write box: */
          fprintf(wr, " %5d %5d  %5d %5d\n", colLO, rowLO, colHI-colLO, rowHI-rowLO);
        }
      fflush(wr);
      free(colmin); free(colmax);
      free(rowmin); free(rowmax);
    }
  }
