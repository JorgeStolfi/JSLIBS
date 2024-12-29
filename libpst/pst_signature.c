/* See pst_signature.h */
/* Last edited on 2024-12-24 19:11:09 by stolfi */ 

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <vec.h>
#include <affirm.h>
#include <float_image.h>
#include <rn.h>

#include <pst_basic.h>
#include <pst_normal_map.h>

#include <pst_signature.h>

#define ltn_normals_debug TRUE

void pst_signature_extract 
  ( image_vec_t *IMGV, /* Scene images under {NF} different light fields. */
    uint32_t maxval,        /* Maxval of original (quantized) images. */                 
    double noise,      /* Additional per-sample noise in images. */                 
    int32_t c,             /* Number of channels in each image. */                      
    int32_t x,             /* Pixel column */                                           
    int32_t y,             /* Pixel row (0 = bottom). */  
    signature_t *sig   /* (OUT) Light signature. */
  )
  { uint32_t NF = IMGV->ne;

    double svar; /* Total variance of sample noise. */
    if (maxval == 0)
      { svar = 0.0; } 
    else
      { double qstep = 1.0/((double)maxval); /* Quantization step. */
        svar = qstep*qstep/12.0;   /* Variance of sample quantization noise. */
      }
    svar += noise*noise;

    /* Extract the samples of channel {c} and compute their Euclidean magnitude {mag}: */
    double magsq = 1.0e-100; /* Channel magnitude squared (plus tiny epsilon to avoid NaNs). */
    double var = 0.0; /* Total variance of image noise in this channel. */
    for (uint32_t i = 0; i < NF; i++)
      { float_image_t *img = IMGV->e[i];
        double p = (double)float_image_get_sample(img, (int32_t)c, x, y); /* Pixel sample value. */
        sig->rin[i] = p; 
        var += svar;
        magsq += p*p;
      }
    double mag = 1.0/sqrt(magsq);
    /* Scale samples of channel {c} by {mag}, adjusting {var} accordingly: */
    for (uint32_t i = 0; i < NF; i++)
      { sig->rin[i] *= mag;
        var *= mag*mag;
      }
    /* Return magnitude and total noise variance of channel: */
    sig->mag = mag;
    sig->var = var;
    if (ltn_normals_debug) 
      { fprintf(stderr, "%3d %5d %5d sig = ", c, x, y);
        pst_signature_print(stderr, NF, sig);
        fprintf(stderr, "\n");
      }
  }
  
void pst_signature_print(FILE *wr, uint32_t NF, signature_t *sig)
  { for (uint32_t i = 0; i < NF; i++) { fprintf(stderr, " %7.4f", sig->rin[i]); }
    fprintf(wr, "  %7.4f", sig->mag);
    fprintf(wr, "  %7.4f", sig->var);
  }

light_table_t *pst_signature_build_table
  ( r2_t *pos,             /* Position in scene images where table is most valid. */
    float_image_t *NRM,    /* Normal map of light field gauge. */ 
    image_vec_t *IMGV,     /* Images of the gauge objects under various light fields. */ 
    bool_t cubify          /* TRUE uses the cube acceleration. */
  )
  {
    uint32_t NF = IMGV->ne;
    demand(NF >= 3, "not enough photos");
    int32_t NC = (int32_t)(IMGV->e[0]->sz[0]);  /* Number of channels in each gauge image. */ 
    int32_t NX = (int32_t)(IMGV->e[0]->sz[1]);  /* Number of columns in each image. */
    int32_t NY = (int32_t)(IMGV->e[0]->sz[2]);  /* Number of rows in each image. */ 
    light_table_t *tab = (light_table_t *)notnull(malloc(sizeof(light_table_t)), "no mem");
    tab->pos = *pos;
    tab->nrm = r3_vec_new(0);  /* {nrm[k]} is the normal of entry {k}. */
    tab->sig = talloc(NC, signature_vec_t);
    for (int32_t c = 0; c < NC; c++) { tab->sig[c] = signature_vec_new(0); }
    uint32_t NE = 0;  /* Number of entries in signature table. */
    /* Scan pixels of light field gauges, gathering the valid table entries. */
    for (int32_t y =0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++) 
          { r3_t nrm;
            bool_t ok = FALSE; /* TRUE iff pixel {(x,y)} is ok, i.e. {nrm} is nonzero. */
            for (uint32_t ax = 0; ax < 3; ax++)
              { float v = float_image_get_sample(NRM, (int32_t)ax, x, y); 
                if (v != 0.0) { ok = TRUE; }
                nrm.c[ax] = (double)v;
              } 
            if (ok)
              { for (int32_t c = 0; c < NC; c++)
                  { signature_vec_expand(&(tab->sig[c]), (vec_index_t)NE);
                    signature_t *sig = &(tab->sig[c].e[NE]);
                    (*sig) = pst_signature_new(NF);
                    pst_signature_extract(IMGV, 0, 0.0, c, x, y, sig);
                  }
                r3_vec_expand(&(tab->nrm), (vec_index_t)NE);
                tab->nrm.e[NE] = nrm;
                NE++;
              }
          }
      }
    demand(NE > 0, "there are no valid pixels in the light gauge images");
    tab->NE = NE;
    /* Trim vectors to exact size: */
    r3_vec_trim(&(tab->nrm), NE);
    for (int32_t c = 0; c < NC; c++) { signature_vec_trim(&(tab->sig[c]), NE); }
    return tab;
  }

void pst_signature_search_table
  ( uint32_t NF,               /* Number of light fields. */
    int32_t NC,               /* Number of color channels. */
    signature_t sig[],    /* Normalized light signature for each channel (size {NC}). */ 
    r2_t pos,             /* Nominal position of {sig} in scene images. */
    uint32_t NG,               /* Number of light-to-normal tables. */
    light_table_t *tab,   /* Light-to-normal tables extracted from gauges (size {NG}). */
    double *dsq,          /* (OUT) Discrepancy squared between {sig} and best match from tables. */
    r3_t *nrm,            /* (OUT) Normal associated to best match in tables. */
    float clr[]           /* (OUT) Intrinisc scene color (size {NC}). */
  )
  {
    uint32_t NS = NF*(uint32_t)NC;
    uint32_t NE = tab->NE;
    
    auto double ediff2(signature_t *sigA, signature_t *sigB);
      /* The expected squared difference between two vectors 
        {a = avgA + delA} and {b = avgB + delB}, of length {NF}; 
        where {avgA = sigA->rin}, {avgB = sigB->rin}, and {delA,delB}
        are independent random vectors with expected mean square
        lengths {sigA->var,sigB->var}, respectively. */
        
    double ediff2(signature_t *sigA, signature_t *sigB)
      { /* We need {E[\abs{(sigA + delA) - (sigB + delB)}^2]}. */
        double sum_d2 = 0.0;
        for (uint32_t j = 0; j < NS; j++) 
          { double dj = sigA->rin[j] - sigB->rin[j];
            sum_d2 += dj*dj;
          }
        return sum_d2 + sigA->var + sigB->var;
      }
      
    /* Brute force search: */
    int32_t kBest = -1;
    double d2Best = +INF;
    { for (int32_t k = 0; k < NE; k++)
        { /* Compute distances for each channel: */
          double sum_wd2 = 0.0;
          double sum_w = 0.0;
          for (int32_t c = 0; c < NC; c++)
            { signature_t *sigAc = &(sig[c]);
              signature_t *sigBck = &(tab->sig[c].e[k]);
              double cd2 = ediff2(sigAc, sigBck);
              double w = 1.0/(sigAc->var + sigBck->var);
              sum_wd2 += w*cd2;
              sum_w += w;
            }
          double d2 = sum_wd2/sum_w;
          if (d2 < d2Best) { kBest = k; d2Best = d2; }
        }
    }

    /* Return normal and mismatch of best match: */
    (*nrm) = tab->nrm.e[kBest];
    (*dsq) = d2Best;

    /* Compute self-color as ratio of per-channel normalizing factors: */
    { for (int32_t c = 0; c < NC; c++)
        { signature_t *sigck = &(tab->sig[c].e[kBest]);
          double cc = sig[c].mag/sigck->mag;
          clr[c] = (float)cc;
        }
    }        
  }

signature_t pst_signature_new(uint32_t NF)
  { signature_t sig;
    sig.rin = rn_alloc(NF);
    sig.mag = 1.0;
    sig.var = 0.0;
    return sig;
  }

void pst_signature_normals_from_photos
  ( uint32_t NG, 
    light_table_t *tab[], 
    image_vec_t *IMGV, 
    uint32_t maxval,
    double noise,
    float_image_t *NRM,
    float_image_t *CLR,
    float_image_t *DIF
  )
  {
    uint32_t NF = IMGV->ne;
    demand(NF >= 3, "too few photos");
    int32_t NC = (int32_t)(IMGV->e[0]->sz[0]);
    int32_t NX = (int32_t)(IMGV->e[0]->sz[1]);
    int32_t NY = (int32_t)(IMGV->e[0]->sz[2]);
    
    /* Check sizes: */
    demand(NRM->sz[0] == 3, "bad num channels in normal map");
    demand(NRM->sz[1] == NX, "bad num columns in normal map");
    demand(NRM->sz[2] == NY, "bad num rows in normal map");
    
    demand(CLR->sz[0] == NC, "bad num channels in color map");
    demand(CLR->sz[1] == NX, "bad num columns in color map");
    demand(CLR->sz[2] == NY, "bad num rows in color map");
    
    demand(DIF->sz[0] == 1, "bad num channels in strangeness map");
    demand(DIF->sz[1] == NX, "bad num columns in strangeness map");
    demand(DIF->sz[2] == NY, "bad num rows in strangeness map");
    
    /* Process scene images, store normal images: */
    r3_t nrm;         /* Computed normal vector of scene surface at some pixel. */
    float clr[NC];    /* Computed intrinsic color of scene surface at some pixel. */
    double dsq;       /* Discrepancy between pixel's signatures and best match from tables. */
    signature_t sig[NC]; /* Light signature at pixel {(x,y)}. */
    for (int32_t c = 0; c < NC; c++) { sig[c] = pst_signature_new(NF); }          
    /* Process scene pixels: */
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { r2_t pos = (r2_t){{ x + 0.5, y + 0.5}}; /* Nominal pixel position. */
            /* Get light signature of scene surface at pixel {(x,y)} in each channel: */
            for (int32_t c = 0; c < NC; c++)
              { pst_signature_extract(IMGV, maxval, noise, c, x, y, &(sig[c])); }
            
            /* Find light field gauge closest to pixel: */
            int32_t gBest = -1;
            { double d2Best = +INF;
              for (int32_t g = 0; g < NG; g++)
                { double d2 = r2_dist_sqr(&pos, &(tab[g]->pos)); 
                  if (d2 < d2Best) { d2Best = d2; gBest = g; }
                }
            }
            /* Map light signature to normal and intrinstic color: */
            pst_signature_search_table(NF, NC, sig, pos, NG, tab[gBest], &dsq, &nrm, clr);
            /* Save discrepancy in {DIF}: */
            float_image_set_sample(DIF, 0, x, y, (float)sqrt(dsq));
            /* Save normal vector in {NRM}: */
            pst_normal_map_set_pixel(NRM, x, y, &nrm);
            /* Save intrinisc color in {CLR}: */
            for (int32_t c = 0; c < NC; c++) { float_image_set_sample(CLR, (int32_t)c, x, y, clr[c]); }
          }
      }
  }

vec_typeimpl(signature_vec_t,signature_vec,signature_t);
