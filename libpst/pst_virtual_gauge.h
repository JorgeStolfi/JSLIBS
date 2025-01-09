#ifndef pst_virtual_gauge_H
#define pst_virtual_gauge_H

/* pst_virtual_gauge.h -- creating a synthetic image of a light gauge. */
/* Last edited on 2025-01-04 04:12:27 by stolfi */

#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <float_image.h>
#include <ellipse_crs.h>
#include <frgb.h>

#include <pst_shading.h>
#include <pst_light.h>
#include <pst_lamp.h>

typedef struct pst_virtual_gauge_data_t
  { ellipse_crs_t E;      /* Ellipse that is the gauge's projection. */
    r3_t view;            /* Viewing direction at the gauge's center. */
    frgb_t albedo;        /* Albedo of gauge in each color channel. */
  } pst_virtual_gauge_data_t;
  /* Describes the geometry of a gauge's projection on an image. */

void pst_virtual_gauge_paint
  ( pst_virtual_gauge_data_t *gd,
    pst_shading_func_t* shading,
    float_image_t *img,
    float_image_t *nrm
  );
  /* Paints into each channel {c} of {img} an image of the virtual
    gauge {gd} with the specified {shading} function.  If {nrm} is not
    {NULL} also paints into it the computed gauge's surface normal.
    
    If a pixel lies entirely outside the gauge's projection, its value
    in {img} is not modeified, and its value in {nrm} is set to
    {(0,0,0)}. If a pixel intersects the gauge's projection, the
    computed gauge color will be mixed into the original {img} pixel
    value width opacity coefficient proportional to the fraction of the
    pixel that is inside the projection. The computed average gauge
    normal will be stored into the {nrm} pixel.
    
    The gauge color and normal in each pixel are obtained by averaging
    the color and normal at a number of sub-sampling points (sampoints)
    inside or near the pixel. The average is taken over the sampoints
    that fall inside the gauge's projection. The normal average is
    adjusted to have unit length.
    
    The image {img} must have 1 to 4 channels. The image of the gauge is
    painted into channels {0..2}. If {img} has 4 channels, the procedure
    writes into channel 3 the opacity of each pixel, thus turning that
    channel a mask image for the projection of the gauge in {img}. */

void pst_virtual_gauge_args_parse(argparser_t *pp, pst_virtual_gauge_data_t *gd);
  /* Parses from the command line the key parameters
    of a virtual gauge on an image, namely the geometty of its outline
    (an ellipse), the viewing direction, and the gauge's albedo. 
    See {pst_virtual_gauge_args_parse_HELP_INFO} below. */
    
#define pst_virtual_gauge_args_parse_HELP \
  ellipse_crs_args_center_HELP " " \
  ellipse_crs_args_radius_HELP " " \
  ellipse_crs_args_stretch_HELP " " \
  pst_virtual_gauge_args_parse_view_HELP \
  pst_virtual_gauge_args_parse_albedo_HELP 
     
#define pst_virtual_gauge_args_parse_HELP_INFO \
  ellipse_crs_args_center_HELP_INFO "\n" \
  "\n" \
  "" ellipse_crs_args_radius_HELP_INFO "\n" \
  "\n" \
  "" ellipse_crs_args_stretch_HELP_INFO "\n" \
  "\n" \
  "" pst_virtual_gauge_args_parse_view_HELP_INFO "\n" \
  "\n" \
  "" pst_virtual_gauge_args_parse_albedo_HELP_INFO 
 
#define pst_virtual_gauge_args_parse_view_HELP \
  "view {IVX} {IVY} {IVZ}"
 
#define pst_virtual_gauge_args_parse_view_INFO \
  "This argument specifies the direction of the observer.  The length" \
  " does not matter, as it is normalized internally to unit length."
   
#define pst_virtual_gauge_args_parse_view_HELP_INFO \
  "  " pst_virtual_gauge_args_parse_view_HELP "\n" \
  "    " pst_virtual_gauge_args_parse_view_INFO ""

#define pst_virtual_gauge_args_parse_albedo_HELP \
  "albedo {IDR} [ {IDG} {IDB} ]" \
 
#define pst_virtual_gauge_args_parse_albedo_INFO \
  "This argument specifies the albedo of the gauge in each color" \
  " channel.  The values {IDR}, {IDG}, {IDB} must be between 0 and 1.  If" \
  " only one number {IDR} is given, it is used for all three channels."
  
#define pst_virtual_gauge_args_parse_albedo_HELP_INFO \
  "  " pst_virtual_gauge_args_parse_albedo_HELP "\n" \
  "    " pst_virtual_gauge_args_parse_albedo_INFO ""

#endif
