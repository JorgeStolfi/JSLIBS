/* See rn_classif_float_image.h. */
/* Last edited on 2024-12-05 10:24:20 by stolfi */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <filefmt.h>
#include <r2.h>
#include <uint16_image.h>
#include <jsrandom.h>
#include <float_image.h>
#include <float_image_paint.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <rn_classif.h>
#include <rn_classif_float_image.h>
        
void rn_classif_float_image_paint_labeler
  ( float_image_t *fim,
    int NA, 
    int NC1, 
    rn_classif_labeler_t *lab1,
    int NC2, 
    rn_classif_labeler_t *lab2, 
    r2_t *ictr, 
    double HV, 
    int NSUB, 
    double sigma,
    frgb_t cmap[]
  )
  {    
    demand(fim->sz[0] == 3, "image is not colored");
    demand(fim->sz[1] == fim->sz[2], "image is not square");
    /* Get image size {NXY} and half-size {HW} of pixel: */
    int NXY = (int)fim->sz[1];
    double HW = HV/NXY;
    
    int ix, iy;
    for (iy = 0; iy < NXY; iy++)
      { for (ix = 0; ix < NXY; ix++) 
          { /* Map {ix,iy} to attribute values in {V} (pixel center): */
            r2_t pctr;
            pctr.c[0] = ictr->c[0] + HV*(2*((double)ix + 0.5)/NXY - 1);
            pctr.c[1] = ictr->c[1] + HV*(2*((double)iy + 0.5)/NXY - 1);
            frgb_t pix = rn_classif_float_image_compute_pixel(NA, NC1, lab1, NC2, lab2, &pctr, HW, NSUB, sigma, cmap);
            float_image_set_pixel(fim, ix, iy, pix.c);
          }
      }
  }
      
void rn_classif_float_image_paint_dataset
  ( float_image_t *fim,
    int NA, 
    int NCD, 
    rn_classif_dataset_t *D,
    double HD[],
    int classD[], 
    r2_t *ictr, 
    double HV,
    bool_t fillH,
    double penH,
    frgb_t cmap[],
    bool_t fillP,
    double radP,
    double penP, 
    int NSUB
  )
  {      
    demand(fim->sz[0] == 3, "image is not colored");
    demand(fim->sz[1] == fim->sz[2], "image is not square");
    /* Get image size {NXY}: */
    int NXY = (int)fim->sz[1];
    
    /* Paint samples: */
    /* fprintf(stderr, "fillH = %c penH = %10.7f\n", "FT"[fillH], penH); */

    int pass;
    for (pass = 0; pass <= 2; pass++)
      { /* Passes: 0 = HD-disk fill, 1 = HD-disk border, 2 = dot fill and border. */ 
        if ((pass == 0) && ((HD == NULL) || (! fillH))) { continue; }
        if ((pass == 1) && ((HD == NULL) || (penH == 0))) { continue; }
        int i, c;
        for (i = 0; i < D->NS; i++)
          { /* Get the attributes {attri} of sample {i} and its class {cli}: */
            double *attri = D->smp[i];
            /* Compute the position of {attri} in the image's domain (pixels from SW corner): */
            double xi = NXY*((attri[0] - ictr->c[0])/(2*HV) + 0.5);
            double yi = NXY*((attri[1] - ictr->c[1])/(2*HV) + 0.5);

            /* Now paint the stuff of the current pass: */
            if (pass < 2) 
              { assert (HD != NULL);
                if (HD[i] != 0) 
                  { /* Show the handicap disk: */
                    double radH = NXY*HD[i]/(2*HV);  /* Disk radius in pixels. */
                    /* fprintf(stderr, "smp[%5d] HD = %10.7f radH = %10.7f\n", i, HD[i], radH); */
                    if (pass == 0)
                      { assert(fillH);
                        /* Paint the disk gray: */
                        double hwdH = 0.0;
                        float vdraw = NAN;
                        frgb_t valH = (frgb_t){{ 0.30f, 0.30f, 0.30f }};
                        bool_t round = TRUE;
                        bool_t diagonal = FALSE;
                        /* Paint the disk with {valH}, no outline: */
                        for (c = 0; c < 3; c++)
                          { float vfill = valH.c[c];
                            (void)float_image_paint_dot
                              ( fim, c, xi, yi, radH, hwdH, 
                                 round, diagonal, vfill, vdraw, NSUB
                              );
                          }
                      }
                    else if (pass == 1)
                      { assert(penH > 0);
                        /* Draw the disk's outline: */
                        double hwdH = penH/2;
                        float vdraw = 0.0;
                        bool_t round = TRUE;
                        bool_t diagonal = FALSE;
                        float vfill = NAN;
                        for (c = 0; c < 3; c++)
                          { (void)float_image_paint_dot
                              ( fim, c, xi, yi, radH, hwdH, 
                                round, diagonal, vfill, vdraw, NSUB
                              );
                          }
                      }
                    else
                      { assert(FALSE); }
                  }
              }
            else if (pass == 2 )
              { /* Paint the dots. */
                double hwdP = penP/2;   /* Pen half-width. */
                float vdraw = 0.0;
                /* Get the point's color {val}: */
                frgb_t val;   /* Nominal dot color. */
                if (! fillP)
                  { /* Paint with clear color: */
                    val = (frgb_t){{ NAN,NAN,NAN }};
                  }
                else if ((classD != NULL) && (cmap != NULL))
                  { /* Color is class color: */
                    int cli = classD[i];
                    assert((cli >= 0) && (cli <= NCD));
                    val = cmap[cli];
                  }
                else
                  { /* Color is black: */
                    val = (frgb_t){{ 0,0,0 }};
                  }
                bool_t round = TRUE;
                bool_t diagonal = FALSE;
                for (c = 0; c < 3; c++)
                  { float vfill = val.c[c];
                    (void)float_image_paint_dot
                      ( fim, c, xi, yi, radP, hwdP, 
                        round, diagonal, vfill, vdraw, NSUB
                      );
                  }
              }
          }
      }
  }

frgb_t rn_classif_float_image_compute_pixel
  ( int NA, 
    int NC1, 
    rn_classif_labeler_t *lab1, 
    int NC2, 
    rn_classif_labeler_t *lab2, 
    r2_t *pctr, 
    double HW, 
    int NSUB, 
    double sigma, 
    frgb_t cmap[]
  )
  {    
    demand(NA == 2, "problem is not two-dimensional");
    int dx, dy, c, t;
    /* Enumerate subsamples in pixel: */
    double sum[3] = { 0,0,0 };
    for (dy = 0; dy < NSUB; dy++)
      { for (dx = 0; dx < NSUB; dx++) 
          { r2_t psmp;
            psmp.c[0] = pctr->c[0] + HW*(2*((double)dx + 0.5)/NSUB - 1);
            psmp.c[1] = pctr->c[1] + HW*(2*((double)dy + 0.5)/NSUB - 1);
            if (sigma > 0)
              { /* Add noise: */
                for (t = 0; t < NA; t++) 
                  { double dit = sigma*fmin(+4.0, fmax(-4.0, dgaussrand()));
                    psmp.c[t] += dit;
                  }
              }
            frgb_t val; /* Sample color. */
            if (lab1 != NULL)
              { /* Classify pixel by labeler 1: */
                int cl1 = lab1(NA, NC1, psmp.c);
                assert((cl1 >= 0) && (cl1 <= NC1));
                val = cmap[cl1];
              }
            else
              { val = (frgb_t){{ 1,1,1 }}; }
            if (lab2 != NULL)
              { /* Classify pixel by labeler 2: */
                int cl2 = lab2(NA, NC2, psmp.c);
                assert((cl2 >= 0) && (cl2 <= NC2));
                /* Use it to modify the saturation and brightness of {val}: */
                frgb_to_HTY(&val);
                double H = val.c[0], T = val.c[1], Y = val.c[2];
                if (cl2 != 0)
                  { /* Reduce the brightness proportionally to {cl2}: */
                    double r = ((double)cl2)/NC2;
                    double Ylim = fmax(0.0, Y - 0.20);
                    Y = Y*(1-r) + Ylim*r;
                  }
                val.c[0] = (float)H; 
                val.c[1] = (float)T; 
                val.c[2] = (float)Y;
                frgb_from_HTY(&val);
              }
            for (c = 0; c < 3; c++) { sum[c] += (double)(val.c[c]); }
          }
      }
    frgb_t res;
    for (c = 0; c < 3; c++) { res.c[c] = (float)(sum[c]/(NSUB*NSUB)); }
    return res;
  }
