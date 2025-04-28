/* See pnmift_image.h */
/* Last edited on 2025-04-24 14:16:39 by stolfi */

#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <jspnm.h>
#include <uint16_image.h>
#include <frgb.h>
#include <affirm.h>

#include <pnmift_image.h>

void pnmift_check_maxval(uint32_t maxval, uint32_t maxmaxval);
  /* Bombs out if {maxval} is zero or greater than {maxmaxval}. */

void pnmift_check_maxval(uint32_t maxval, uint32_t maxmaxval)
  { if ((maxval < 1) || (maxval > maxmaxval)) 
      { fprintf(stderr, "maxval = %d\n", maxval); 
        demand(FALSE, "bad maxval");
      }
  }

void pnmift_set_values_from_image(uint16_image_t *img, ift_graph_t *G, frgb_t rgb[])
  { double fmaxval = (double)(img->maxval);
    demand(img->cols == G->cols, "cols mismatch"); 
    demand(img->rows == G->rows, "rows mismatch");
    demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    /* Initialize unused channels: */
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col*(int32_t)img->chns]);
          frgb_t yy;
          for (uint32_t chn = 0; chn < 3; chn++) 
            { if (chn < img->chns)
                { yy.c[chn] = (float)(((double)(*pm))/fmaxval); pm++; }
              else
                { yy.c[chn] = yy.c[0]; }
            }
          int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          assert(pg->col == col);
          assert(pg->row == row);
          rgb[ig] = yy;
        }
  }

void pnmift_set_labels_from_image(uint16_image_t *seed_img, ift_graph_t *G, uint16_t label[], bool_t verbose)
  { demand(seed_img->chns == 1, "seed image must be monochromatic");
    demand(seed_img->cols == G->cols, "cols mismatch");
    demand(seed_img->rows == G->rows, "rows mismatch");
    demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { uint32_t pm = seed_img->smp[row][col];
          demand(pm <= pnmift_MAX_LABEL, "seed label too big");
          int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          assert(pg->col == col);
          assert(pg->row == row);
          label[ig] = (uint16_t)pm;
          if (verbose && (pm != pnmift_NON_SEED_LABEL))
            { fprintf(stderr, "  ( %5d %5d ) = %d\n", (ift_pixel_index_t)col, (ift_pixel_index_t)row, pm); }
        }
  }

uint16_image_t *pnmift_get_cost_image(ift_graph_t *G, ift_path_cost_t maxcost, uint16_t maxval)
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    uint16_image_t *img = uint16_image_new((uint32_t)G->cols, (uint32_t)G->rows, 1);
    pnmift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);
    img->maxval = maxval;
    double scale = (maxcost <= 0 ? 1.0 : ((double)(maxval - 1))/maxcost);
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col]);
          int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          if (pg->C == +INFINITY)
            { (*pm) = img->maxval; }
          else
            { double y = scale*pg->C;
              if (y < 0) 
                { (*pm) = 0; }
              else if (y > (double)(img->maxval-1)) 
                { (*pm) = (uint16_t)(img->maxval-1); }
              else
                { (*pm) = (uint16_t)(y + 0.5); }
            }
        }
    return img;
  }

uint16_image_t *pnmift_get_pred_image(ift_graph_t *G)
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    uint16_image_t *img = uint16_image_new(G->cols, G->rows, 1);
    img->maxval = PNM_FILE_MAX_MAXVAL;
    demand(img->maxval >= 65535, "maxval too small"); 
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col]);
          int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          ift_node_t *pp = (ift_node_t *)pg->P;
          int32_t dcol, drow;
          if (pp == NULL)
            { dcol = 0; drow = 0; }
          else
            { dcol = (int32_t)pp->col - col; drow = (int32_t)pp->row - row; }
          demand((dcol >= -127) && (dcol <= 127), "pred dcol out of range"); 
          demand((drow >= -127) && (drow <= 127), "pred drow out of range"); 
          (*pm) = (uint16_t)(((drow + 128) << 8) | (dcol + 128));
        }
    return img;
  }

uint16_image_t *pnmift_get_label_image(ift_graph_t *G, uint16_t label[], uint16_t maxval)
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    uint16_image_t *img = uint16_image_new(G->cols, G->rows, 1);
    pnmift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);
    img->maxval = maxval;
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col]);
          int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          ift_node_t *pr = pg->R; /* Root of forest containing {pg}. */
          demand(pr != NULL, "null root"); 
          int32_t ir = ift_node_index(G, pr->col, pr->row);
          uint16_t lab = label[ir];
          /* fprintf(stderr, " ( %5d %5d ) = %5d --> ( %5d %5d ) = %5d\n", col,row,label[ig],pr->col,pr->row,label[ir]); */
          demand(lab <= maxval, "label too big"); 
          (*pm) = lab;
        }
    return img;
  }

uint16_image_t *pnmift_get_root_image(ift_graph_t *G)
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    uint16_image_t *img = uint16_image_new(G->cols, G->rows, 1);
    img->maxval = 255;
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col]);
          int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          demand((pg->P == NULL) == (pg->R == pg), "P/R inconsistency"); 
          (*pm) = (pg->P == NULL ? img->maxval : 0);
        }
    return img;
  }

uint16_image_t *pnmift_get_spread_image(ift_graph_t *G, frgb_t rgb[], uint32_t chns, uint16_t maxval)
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    uint16_image_t *img;
    double fmaxval = (double)maxval;
    /* Choose image type: */
    pnmift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);

    /* Allocate image and fill it out: */
    img = uint16_image_new(G->cols, G->rows, chns);
    img->maxval = maxval;
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          ift_node_t *pr = pg->R;
          demand(pr != NULL, "R is null"); 
          demand((pr == pg) == (pg->P == NULL), "P/R inconsistency"); 
          int32_t ir = ift_node_index(G, pr->col, pr->row);
          frgb_t yy = rgb[ir];
          uint16_t *pm = &(img->smp[row][col*(int32_t)img->chns]);
          for (uint32_t chn = 0; chn < chns; chn++)
            { double fv = yy.c[chn];
              int32_t pv;
              if (fv <= 0.0) 
                { pv = 0; }
              else if (fv >= 1.0)
                { pv = maxval; }
              else
                { pv = (int32_t)floor(0.5 + fv*fmaxval); }
              (*pm) = (uint16_t)pv;
              pm++;
            }
        }
    return img;
  }

uint16_image_t *pnmift_get_single_label_image
  ( ift_graph_t *G, 
    uint16_t label[], 
    uint16_t lab, 
    frgb_t rgb[],
    uint16_t bg[], 
    uint32_t chns,
    uint16_t maxval
  )
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    double fmaxval = (double)maxval;
    /* Choose imag etype: */
    pnmift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);

    /* Allocate image: */
    uint16_image_t *img = uint16_image_new(G->cols, G->rows, chns);
    img->maxval = maxval;

    /* Fill image: */
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          frgb_t yy = rgb[ig];
          demand(pg->R != NULL, "null root"); 
          ift_node_t *pr = pg->R;
          ift_node_index_t ir = ift_node_index(G, pr->col, pr->row);
          uint16_t prlab = label[ir];
          /* Set the {chns} channels of {img.smp[row][col]} from {yy} or {bg}: */
          uint16_t *pm = &(img->smp[row][col*(int32_t)img->chns]);
          for (uint32_t chn = 0; chn < chns; chn++)
            { int32_t pv;
              if (prlab == lab) 
                { double fv = yy.c[chn];
                  if (fv <= 0.0) 
                    { pv = 0; }
                  else if (fv >= 1.0)
                    { pv = maxval; }
                  else
                    { pv = (int32_t)floor(0.5 + fv*fmaxval); }
                }
              else
                { pv = bg[chn]; }
              (*pm) = (uint16_t)pv;
	      pm++;
            }
        }
    return img;
  }

void pnmift_write_boxes(FILE *wr, ift_graph_t *G, uint16_t label[], uint16_t maxval, int32_t margin)
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    
    /* Max and min pixel indices for each label: */
    int32_t *colmin = (int32_t *)notnull(malloc((maxval+1)*sizeof(int32_t)), "out of mem");
    int32_t *colmax = (int32_t *)notnull(malloc((maxval+1)*sizeof(int32_t)), "out of mem");
    int32_t *rowmin = (int32_t *)notnull(malloc((maxval+1)*sizeof(int32_t)), "out of mem");
    int32_t *rowmax = (int32_t *)notnull(malloc((maxval+1)*sizeof(int32_t)), "out of mem");
    
    pnmift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);

    /* Initialize max/min accumulators: */
    { for (uint32_t lab = 0; lab <= maxval; lab++) 
        { colmin[lab] = rowmin[lab] = 1000000; 
          colmax[lab] = rowmax[lab] = 0;
        }
    }
      
    /* Scan image and collect max/min indices per label: */
    for (ift_pixel_index_t row = 0; row < G->rows; row++)
      for (ift_pixel_index_t col = 0; col < G->cols; col++)
        { int32_t ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          demand(! (pg->R == NULL), "null root"); 
          ift_node_t *pr = pg->R;
          ift_node_index_t ir = ift_node_index(G, pr->col, pr->row);
          uint16_t lab = label[ir];
          demand(lab <= maxval, "label too big"); 
          if (col < colmin[lab]) { colmin[lab] = col; }
          if (col > colmax[lab]) { colmax[lab] = col; }
          if (row < rowmin[lab]) { rowmin[lab] = row; }
          if (row > rowmax[lab]) { rowmax[lab] = row; }
        }
      
    /* Print region boxes: */
    { for (uint32_t lab = 0; lab <= maxval; lab++) 
        { int32_t colLO = colmin[lab], colHI = colmax[lab] + 1;
          int32_t rowLO = rowmin[lab], rowHI = rowmax[lab] + 1;
          fprintf(wr, "%5d ", lab);
          /* Apply margin expansion/contraction: */
          if ((colLO < colHI) && (rowLO < rowHI))
            { colLO -= margin; colHI += margin;
              rowLO -= margin; rowHI += margin;
              if (colLO < 0) { colLO = 0; }
              if (colHI > G->cols) { colHI = (int32_t)G->cols; }
              if (rowLO < 0) { rowLO = 0; }
              if (rowHI > G->rows) { rowHI = (int32_t)G->rows; }
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
