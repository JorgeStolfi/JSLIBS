/* See {uint16_image_quantize_floyd.h}.  */
/* Last edited on 2024-12-26 19:50:24 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <jspnm.h>
#include <jsrandom.h>
#include <argparser.h>
#include <affirm.h>
#include <bool.h>
#include <uint16_image.h>
#include <uint16_image_match.h>
#include <uint16_color_table.h>
#include <uint16_color_tree.h>

#include <uint16_image_quantize_floyd.h>

void uint16_image_quantize_floyd
  ( uint16_image_t *img, 
    uint16_color_table_t *chv, 
    uint16_t mapmaxval,
    uint32_t maxColors,
    uint16_image_match_proc_t *match
  )
  { int32_t chns = img->chns; assert(chns == 3);
    int32_t cols = img->cols;
    int32_t rows = img->rows;
    assert(mapmaxval <= PNM_FILE_MAX_MAXVAL);
    
#define FS_SCALE 1024
#define FS_WTA 7
#define FS_WTB 3
#define FS_WTC 5
#define FS_WTD 1
    int32_t* thisrerr = talloc(img->cols + 2, int32_t);
    int32_t* nextrerr = talloc(img->cols + 2, int32_t);
    int32_t* thisgerr = talloc(img->cols + 2, int32_t);
    int32_t* nextgerr = talloc(img->cols + 2, int32_t);
    int32_t* thisberr = talloc(img->cols + 2, int32_t);
    int32_t* nextberr = talloc(img->cols + 2, int32_t);
    int32_t* temperr;
    int32_t row, col;
    uint16_t maxval = img->maxval;
    int32_t maxcor = 2*(int32_t)maxval;
    int32_t mincor = -1*(int32_t)maxval;
    int32_t fs_direction = 1;
    bool_t addtohash = 1;
    uint16_image_RGB_table_node_t **cht = uint16_image_RGB_table_alloc();
    int32_t sr, sg, sb, err;
    uint32_t new;
    float coef = FS_SCALE*((float)maxval)/((float)mapmaxval);

    /* Initialize Floyd-Steinberg error vectors with randoms in [-1..+1]. */
    srand((int32_t) 46157);
    for (int32_t col = 0; col < cols + 2; col++)
      { thisrerr[col] = int32_abrandom(-FS_SCALE, +FS_SCALE-1);
        thisgerr[col] = int32_abrandom(-FS_SCALE, +FS_SCALE-1);
        thisberr[col] = int32_abrandom(-FS_SCALE, +FS_SCALE-1);
      }

    /* Scan the image rows: */
    for (row = 0; row < rows; row++)
      { /* Clear the {nextXerr} arrays: */
        for (col = 0; col < cols + 2; col++)
          { nextrerr[col] = nextgerr[col] = nextberr[col] = 0; }
        /* Define the initial and final column */
        int32_t limitcol;
        if (fs_direction)
          { col = 0; limitcol = cols; }
        else
          { col = cols - 1; limitcol = -1; }
        /* Scan columns and propagate errors */
        uint16_t *pp = &(img->smp[row][col*chns]);
        do
          { /* Use Floyd-Steinberg errors to adjust actual color. */
            sr = ((int32_t)pp[0] + thisrerr[col + 1]) / FS_SCALE;
            sg = ((int32_t)pp[1] + thisgerr[col + 1]) / FS_SCALE;
            sb = ((int32_t)pp[2] + thisberr[col + 1]) / FS_SCALE;
            if (sr < mincor) { sr = mincor; } if (sr > maxcor) { sr = maxcor; }
            if (sg < mincor) { sg = mincor; } if (sg > maxcor) { sg = maxcor; }
            if (sb < mincor) { sb = mincor; } if (sb > maxcor) { sb = maxcor; }
            /* Choose replacement color* */
            ppm_pixel_t *qq = uint16_image_RGB_table_choose_color
              ( sr, sg, sb, maxval, 
                chv, newcolors, mapmaxval,
                match, cht, &addtohash
              );
            /* Replace pixel in image: */
            pp[0] = qq->c[0];
            pp[1] = qq->c[1];
            pp[2] = qq->c[2];
            /* Propagate Floyd-Steinberg error terms. */
            if (fs_direction)
              { new = (uint32_t)(coef*qq->c[0]+0.5);
                err = (sr * FS_SCALE - (int32_t)new);
                thisrerr[col + 2] += (err * FS_WTA) / 16;
                nextrerr[col    ] += (err * FS_WTB) / 16;
                nextrerr[col + 1] += (err * FS_WTC) / 16;
                nextrerr[col + 2] += (err * FS_WTD) / 16;
                new = (uint32_t)(coef*qq->c[1]+0.5);
                err = (sg * FS_SCALE - (int32_t)new);
                thisgerr[col + 2] += (err * FS_WTA) / 16;
                nextgerr[col    ] += (err * FS_WTB) / 16;
                nextgerr[col + 1] += (err * FS_WTC) / 16;
                nextgerr[col + 2] += (err * FS_WTD) / 16;
                new = (uint32_t)(coef*qq->c[2]+0.5);
                err = (sb * FS_SCALE - (int32_t)new);
                thisberr[col + 2] += (err * FS_WTA) / 16;
                nextberr[col    ] += (err * FS_WTB) / 16;
                nextberr[col + 1] += (err * FS_WTC) / 16;
                nextberr[col + 2] += (err * FS_WTD) / 16;
                col++; pp += chns;
              }
            else
              { new = (uint32_t)(coef*qq->c[0]+0.5);
                err = (sr * FS_SCALE - (int32_t)new);
                thisrerr[col    ] += (err * FS_WTA) / 16;
                nextrerr[col + 2] += (err * FS_WTB) / 16;
                nextrerr[col + 1] += (err * FS_WTC) / 16;
                nextrerr[col    ] += (err * FS_WTD) / 16;
                new = (uint32_t)(coef*qq->c[1]+0.5);
                err = (sg * FS_SCALE - (int32_t)new);
                thisgerr[col    ] += (err * FS_WTA) / 16;
                nextgerr[col + 2] += (err * FS_WTB) / 16;
                nextgerr[col + 1] += (err * FS_WTC) / 16;
                nextgerr[col    ] += (err * FS_WTD) / 16;
                new = (uint32_t)(coef*qq->c[2]+0.5);
                err = (sb * FS_SCALE - (int32_t)new);
                thisberr[col    ] += (err * FS_WTA) / 16;
                nextberr[col + 2] += (err * FS_WTB) / 16;
                nextberr[col + 1] += (err * FS_WTC) / 16;
                nextberr[col    ] += (err * FS_WTD) / 16;
                col--; pp -= chns;                
              }
          }
        while (col != limitcol);

        temperr = thisrerr;
        thisrerr = nextrerr;
        nextrerr = temperr;
        temperr = thisgerr;
        thisgerr = nextgerr;
        nextgerr = temperr;
        temperr = thisberr;
        thisberr = nextberr;
        nextberr = temperr;
        fs_direction = ! fs_direction;
      }
    img->maxval = (uint16_t)mapmaxval;
  }
