/* See uint16_image_Canon_EOS50D.h */
/* Last edited on 2017-06-20 20:51:43 by stolfilocal */

#include <stdlib.h>
#include <assert.h>

#include <affirm.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_Pike_F100.h>

#include <uint16_image_Canon_EOS50D.h>

uint16_image_t *pnm_Canon_EOS50D_raw16_debayer(uint16_image_t *img, bool_t squeeze, bool_t verbose)  
  {
    demand(img->chns == 1, "input must be single-channel");
    uint16_image_t *omg;
    if (! squeeze)
      { /* Retain the original size, just separate the channels: */
        int rows = img->rows;
        int cols = img->cols;
        omg = uint16_image_new(cols, rows, 3);
        omg->maxval = img->maxval;
        int col, row;
        for (row = 0; row < rows; row++)
          { for (col = 0; col < cols; col++)
              { /* Decide which channel {chn} this sample belongs to: */
                int chn = uint16_image_Pike_F100_sample_channel(col, row);
                /* Copy sample, set other channels to 0: */
                uint16_t smp = img->smp[row][col];
                uint16_t *po = &(omg->smp[row][3*col]);
                int k;
                for (k = 0; k < 3; k++) 
                  { (*po) = (uint16_t)(k == chn ? smp : 0); po++; }
              }
          }
      }
    else
      { /* Allocate a half-size 3-channel image: */
        omg = uint16_image_new(img->cols/2, img->rows/2, 3);
        omg->maxval = img->maxval;
        /* Fill it with the merged channels: */
        int col,row;
        for (row = 0; row < omg->rows; row++)
          { 
            uint16_t *R0p = &(img->smp[2*row][0]);
            uint16_t *G0p = &(img->smp[2*row][1]);
            uint16_t *G1p = &(img->smp[2*row+1][0]);
            uint16_t *B0p = &(img->smp[2*row+1][1]);
	    
	    uint16_t *op = omg->smp[row];
            for (col = 0; col < omg->cols; col++)
              { (*op) = (*R0p); op++; R0p++; R0p++;
                int G0v = (int)(*G0p); G0p++; G0p++;
                int G1v = (int)(*G1p); G1p++; G1p++;
                /* Apply a simple dither matrix to avoid downwards bias in rounding: */
                int dither = (col + row) % 2;
                int Gv = (G0v + G1v + dither)/2;
                assert((Gv >= 0) && (Gv <= omg->maxval));
                (*op) = (uint16_t)Gv; op++;
                (*op) = (*B0p); op++; B0p++; B0p++;
              }
          }
      }
    return omg;
  }  
               
