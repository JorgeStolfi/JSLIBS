/* See {multifok_frame.h}. */
/* Last edited on 2025-04-11 14:50:19 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_map_channels.h>

#include <multifok_image.h>
#include <multifok_image_read.h>
#include <multifok_image_write.h>
#include <multifok_weight.h>

#include <multifok_frame.h>
  
#define DASHES "----------------------------------------------------------------------"
     
multifok_frame_t *multifok_frame_from_images
  ( float_image_t *sVal,    /* Simulated camera image. */
    float_image_t *shrp,    /* Actual sharpness map. */
    float_image_t *hAvg,    /* Actual scene {Z} average map. */
    float_image_t *hDev,    /* Actual scene {Z} deviation map. */
    float_image_t *sNrm,    /* Surface normal direction. */
    double zFoc,            /* Nominal {Z} of focus plane in each image. */
    double zDep             /* Nominal depth of focus. */
  )
  {
    multifok_frame_t *frame = talloc(1, multifok_frame_t);
 
    frame->sVal = sVal;
    frame->shrp = shrp; 
    frame->hAvg = hAvg; 
    frame->hDev = hDev; 
    frame->sNrm = sNrm; 

    float_image_get_size(sVal, &(frame->NC), &(frame->NX), &(frame->NY));
    float_image_check_size(frame->shrp, 1, frame->NX, frame->NY, "bad sharp image");
    float_image_check_size(frame->hAvg, 1, frame->NX, frame->NY, "bad height map image");
    float_image_check_size(frame->hDev, 1, frame->NX, frame->NY, "bad height deviation map image");
    float_image_check_size(frame->sNrm, 3, frame->NX, frame->NY, "bad normal map image");
    
    frame->zFoc = zFoc;
    frame->zDep = zDep;
    
    return frame;
  }
  
void multifok_frame_free(multifok_frame_t *frame)
  {
    float_image_free(frame->sVal);
    float_image_free(frame->hAvg);
    float_image_free(frame->hDev);
    float_image_free(frame->shrp);
    float_image_free(frame->sNrm);
    free(frame);
  }

multifok_frame_t *multifok_frame_read
  ( char *frameFolder,
    bool_t gray,
    double zFoc,
    double zDep,
    double hMin,
    double hMax
  )
  {
    float_image_t *sVal_in = multifok_image_read_scene_view(frameFolder); 
    int32_t NC, NX, NY;
    float_image_get_size(sVal_in, &(NC), &(NX), &(NY));
    float_image_t *sVal = NULL;
    if (gray & (NC != 1)) 
      { sVal = float_image_new(1, NX, NY);
        /* Discard {U,V} channels. */
        float_image_map_channels_RGB_to_YUV(sVal_in, sVal); 
        float_image_free(sVal_in);
        NC = 1;
      }
    else
      { sVal = sVal_in; }
    float_image_t *hAvg = multifok_image_read_height_average(frameFolder, hMin, hMax); 
    float_image_t *hDev = multifok_image_read_height_deviation(frameFolder, (hMax-hMin)/2); 
    float_image_t *sNrm = multifok_image_read_normal_average(frameFolder); 
    float_image_t *shrp = multifok_image_read_sharpness(frameFolder); 
    
    multifok_frame_t *frame = multifok_frame_from_images
      ( sVal, shrp, hAvg, hDev, sNrm, zFoc, zDep );
    
    return frame;
  }
 
void multifok_frame_write
  ( multifok_frame_t *frame,
    char *frameFolder,
    double hMin,
    double hMax
  )
  {
    mkdir(frameFolder, 0755); /* "-rwxr-xr-x" */
    multifok_image_write_scene_view(frame->sVal, frameFolder);
    multifok_image_write_height_average(frame->hAvg, frameFolder, hMin, hMax);
    multifok_image_write_height_deviation(frame->hDev, frameFolder, hMax-hMin);
    multifok_image_write_normal_average(frame->sNrm, frameFolder);
    multifok_image_write_sharpness(frame->shrp, frameFolder);
    if (frame->zDep == +INF)
      { /* Sharp frame. Write the FNI versions of the height and normal maps: */
      
        /* Compute a weight map: */
        double dhMax = 0.05*(hMax-hMin); /* Max height deviation in a pixel. */
        double dnMax = 0.1; /* Max normal difference between adh pixels. */
        float_image_t *sWht = multifok_weight_from_height_and_normal(frame->hDev, dhMax, frame->sNrm, dnMax);

        /* Add weight to height map and write: */
        float_image_t *hAvgWht = multifok_image_set_weight_channel(frame->hAvg, sWht, 1, TRUE);
        multifok_image_write_fni_height_average(hAvgWht, frameFolder);

        /* Add weight to normal map and write: */
        float_image_t *sNrmWht = multifok_image_set_weight_channel(frame->sNrm, sWht, 3, TRUE);
        multifok_image_write_fni_normal_average(sNrmWht, frameFolder);

        float_image_free(sWht);
        float_image_free(hAvgWht);
        float_image_free(sNrmWht);
      }
  }
  
#define multifok_frame_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

