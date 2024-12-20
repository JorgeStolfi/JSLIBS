/* See pnmift_image.h */
/* Last edited on 2024-12-05 10:29:12 by stolfi */

#include <uint16_image.h>
#include <math.h>
#include <assert.h>

#include <jspnm.h>
#include <uint16_image.h>
#include <frgb.h>
#include <affirm.h>

#include <pnmift_image.h>

void pnmift_check_maxval(unsigned int maxval, unsigned int maxmaxval);
  /* Bombs out if {maxval} is zero or greater than {maxmaxval}. */

void pnmift_check_maxval(unsigned int maxval, unsigned int maxmaxval)
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
    int col, row, chn; 
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col*img->chns]);
          frgb_t yy;
          for (chn = 0; chn < 3; chn++) 
            { if (chn < img->chns)
                { yy.c[chn] = (float)(((double)(*pm))/fmaxval); pm++; }
              else
                { yy.c[chn] = yy.c[0]; }
            }
          int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
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
    int col, row;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { unsigned int pm = seed_img->smp[row][col];
          demand(pm <= pnmift_MAX_LABEL, "seed label too big");
          int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
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
    uint16_image_t *img = uint16_image_new(G->cols, G->rows, 1);
    pnmift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);
    img->maxval = maxval;
    double scale = (maxcost <= 0 ? 1.0 : ((double)(maxval - 1))/maxcost);
    int col, row;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col]);
          int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
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
    int col, row, dcol, drow;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col]);
          int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          ift_node_t *pp = (ift_node_t *)pg->P;
          if (pp == NULL)
            { dcol = 0; drow = 0; }
          else
            { dcol = pp->col - col; drow = pp->row - row; }
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
    int col, row;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col]);
          int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          ift_node_t *pr = pg->R; /* Root of forest containing {pg}. */
          demand(pr != NULL, "null root"); 
          int ir = ift_node_index(G, pr->col, pr->row);
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
    int col, row;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { uint16_t *pm = &(img->smp[row][col]);
          int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          demand((pg->P == NULL) == (pg->R == pg), "P/R inconsistency"); 
          (*pm) = (pg->P == NULL ? img->maxval : 0);
        }
    return img;
  }

uint16_image_t *pnmift_get_spread_image(ift_graph_t *G, frgb_t rgb[], int chns, uint16_t maxval)
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    uint16_image_t *img;
    int col, row, chn;
    double fmaxval = (double)maxval;
    /* Choose image type: */
    pnmift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);

    /* Allocate image and fill it out: */
    img = uint16_image_new(G->cols, G->rows, chns);
    img->maxval = maxval;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          ift_node_t *pr = pg->R;
          demand(pr != NULL, "R is null"); 
          demand((pr == pg) == (pg->P == NULL), "P/R inconsistency"); 
          int ir = ift_node_index(G, pr->col, pr->row);
          frgb_t yy = rgb[ir];
          uint16_t *pm = &(img->smp[row][col*img->chns]);
          for (chn = 0; chn < chns; chn++)
            { double fv = yy.c[chn];
              int pv;
              if (fv <= 0.0) 
                { pv = 0; }
              else if (fv >= 1.0)
                { pv = maxval; }
              else
                { pv = (int)floor(0.5 + fv*fmaxval); }
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
    int chns,
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
    int col, row, chn;
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
          ift_node_t *pg = &(G->node[ig]);
          frgb_t yy = rgb[ig];
          demand(pg->R != NULL, "null root"); 
          ift_node_t *pr = pg->R;
          ift_node_index_t ir = ift_node_index(G, pr->col, pr->row);
          uint16_t prlab = label[ir];
          /* Set the {chns} channels of {img.smp[row][col]} from {yy} or {bg}: */
          uint16_t *pm = &(img->smp[row][col*img->chns]);
          for (chn = 0; chn < chns; chn++)
            { int pv;
              if (prlab == lab) 
                { double fv = yy.c[chn];
                  if (fv <= 0.0) 
                    { pv = 0; }
                  else if (fv >= 1.0)
                    { pv = maxval; }
                  else
                    { pv = (int)floor(0.5 + fv*fmaxval); }
                }
              else
                { pv = bg[chn]; }
              (*pm) = (uint16_t)pv;
	      pm++;
            }
        }
    return img;
  }

void pnmift_write_boxes(FILE *wr, ift_graph_t *G, uint16_t label[], uint16_t maxval, int margin)
  { demand((G->rows <= pnmift_MAX_ROWS) && (G->cols <= pnmift_MAX_COLS), "image is too big");
    int col, row;
    
    /* Max and min pixel indices for each label: */
    int *colmin = (int *)notnull(malloc((maxval+1)*sizeof(int)), "out of mem");
    int *colmax = (int *)notnull(malloc((maxval+1)*sizeof(int)), "out of mem");
    int *rowmin = (int *)notnull(malloc((maxval+1)*sizeof(int)), "out of mem");
    int *rowmax = (int *)notnull(malloc((maxval+1)*sizeof(int)), "out of mem");
    
    pnmift_check_maxval(maxval, PNM_FILE_MAX_MAXVAL);

    /* Initialize max/min accumulators: */
    { int lab;
      for (lab = 0; lab <= maxval; lab++) 
        { colmin[lab] = rowmin[lab] = 1000000; 
          colmax[lab] = rowmax[lab] = 0;
        }
    }
      
    /* Scan image and collect max/min indices per label: */
    for (row = 0; row < G->rows; row++)
      for (col = 0; col < G->cols; col++)
        { int ig = ift_node_index(G, (ift_pixel_index_t)col, (ift_pixel_index_t)row);
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
