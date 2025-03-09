/* See pst_height_map.h */
/* Last edited on 2025-03-04 18:28:48 by stolfi */

#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsprintf.h>
#include <float_image.h>
#include <pst_map.h>
#include <pst_vertex_map_shrink.h>
#include <pst_vertex_map_expand.h>
#include <pst_map_compare.h>

#include <pst_height_map.h>

float_image_t *pst_height_map_shrink(float_image_t *IZ, double scale)
  {
    bool_t debug = FALSE;
    
    int32_t NC, NXI, NYI;
    float_image_get_size(IZ, &NC, &NXI, &NYI);
    demand((NC == 1) || (NC == 2), "height map must have 1 or 2 channels");

    int32_t NXJ = NXI/2+1;
    int32_t NYJ = NYI/2+1;
    
    /* Create the output image {JZ}: */
    if (debug) { fprintf(stderr, "creating shrunk image %d×%d ... \n", NXJ, NYJ); }
    float_image_t *JZ = pst_vertex_map_shrink(IZ,1,NXJ,NYJ,scale);
    return JZ;
  }

float_image_t *pst_height_map_expand(float_image_t *JZ, int32_t NXI, int32_t NYI, double scale)
  { 
    bool_t debug = FALSE;
  
    int32_t NC, NXJ, NYJ;
    float_image_get_size(JZ, &NC, &NXJ, &NYJ);
    demand((NC == 1) || (NC == 2), "input height map must have 1 or 2 channels");
    
    demand(NXJ == NXI/2 + 1, "inconsistent {X} sizes");
    demand(NYJ == NYI/2 + 1, "inconsistent {Y} sizes");
    
    if (debug) { fprintf(stderr, "creating expanded image %d×%d ... \n", NXJ, NYJ); }
    float_image_t *IZ = pst_vertex_map_expand(JZ,1,NXI,NYI,scale);
    return IZ;
  }

void pst_height_map_perturb(float_image_t *A, int32_t wch, double relNoise, double absNoise)
  { bool_t debug = TRUE;
    demand(isfinite(relNoise) && (relNoise >= 0), "invalid {relNoise}");
    demand(isfinite(absNoise) && (absNoise >= 0), "invalid {absNoise}");
    
    int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);

    for (int32_t c = 0; c < NC; c++)
      { if (c != wch)
          { /* Get original sample range in channel {c}: */
            float vmin = +INF, vmax = -INF;
            float_image_update_sample_range(A, c, &vmin, &vmax);
            if (debug) { fprintf(stderr, "  channel %d original range [%.6f _ %.6f]\n", c, vmin, vmax); }
            if (vmin > vmax) { vmin = vmax = 0; }

            /* Decide magnitude {absNoise} of perturbation: */
            double mag = hypot(absNoise, relNoise*(vmax - vmin));
            if (debug) { fprintf(stderr, "  perturbing channel %d with [%.6f _ %.6f]\n", c, -mag, +mag); }
            pst_height_map_perturb_channel(A, c, mag);
          }
      }
    pst_map_ensure_pixel_consistency(A, wch);
  }

void pst_height_map_perturb_channel(float_image_t *A, int32_t c, double mag)
  {
    int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);

    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float *vxy = float_image_get_sample_address(A, c, x, y);
            if (isfinite(*vxy))
              { /* Generate a sum {pert} of threelow-freq  waves: */
                double t1 = (+2.0*x)/NX + (+1.0*y)/NY + 0.25*(1+c); double v1 = sin(2*M_PI*t1);
                double t2 = (-1.0*x)/NX + (-3.0*y)/NY + 0.50*(1+c); double v2 = sin(2*M_PI*t2);
                double t3 = (-2.0*x)/NX + (+2.0*y)/NY + 0.75*(1+c); double v3 = sin(2*M_PI*t3);
                double pert = v1 + v2 + v3;

                /* Squash the perturbation {pert} from {[-3 _ +3]} to {[-1 _ +1]}: */
                pert = sin(0.5*M_PI*pert/3);

                /* Modify the sample: */
                float pvxy = (float)((*vxy) + mag*pert);
                if (isfinite(pvxy)) { (*vxy) = pvxy; }
              }
          }
      }
  }

void pst_height_map_analyze_and_write
  ( char *filePrefix,
    int32_t level,
    int32_t iter,
    double change,
    float_image_t *Z,
    float_image_t *R, 
    bool_t verbose
  )
  {
    bool_t debug = TRUE;
    int32_t indent = (level <= -1 ? 0 : 2*level+2);
    if (debug) { fprintf(stderr, "%*sentering {%s} iter = %d\n", indent, "", __FUNCTION__, iter); }
    
    int32_t NC_Z, NX, NY;
    float_image_get_size(Z, &NC_Z, &NX, &NY);
    demand((NC_Z == 1) || (NC_Z == 2), "input height map must have 1 or 2 channels");
    
    if (R != NULL) 
      { float_image_check_size(R, -1, NX, NY, "height map size mismatch");
        int32_t NC_R = (int32_t)R->sz[0];
        demand((NC_R == 1) || (NC_R == 2), "reference height map must have 1 or 2 channels");
      }
    
    char *tagA = "Z";
    char *tagB = NULL;
    char *tagE = (R != NULL ? "E" : NULL);
    
    double changev[2]; changev[0] = change; changev[1] = 0;
    int32_t NC = 2;
    pst_map_compare_analyze_and_write(Z, R, NC, changev, filePrefix, level, iter, tagA, tagB, tagE, verbose); 
    if (debug) { fprintf(stderr, "%*sexiting {%s}\n", indent, "", __FUNCTION__); }
  }
