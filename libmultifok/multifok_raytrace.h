/* Generic raytracing for the {multifok} library. */
/* Last edited on 2025-02-05 14:59:47 by stolfi */

#ifndef multifok_raytrace_H
#define multifok_raytrace_H

#include <stdint.h>

#include <bool.h>
#include <i2.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>

#include <multifok_frame.h>
#include <multifok_sampling.h>

#define ix_DEBUG 10
#define iy_DEBUG 17
  /* Default pixel indices for which debugging info is to be printed out. */

typedef void multifok_raytrace_proc_t
  ( r3_t *pRay,
    r3_t *dRay,
    bool_t debug,
    r3_t *pHit_P,
    r3_t *sNrm_P,
    int32_t NC,
    float sVal[]
  );
  /* Type of a procedure suitable for the {trace_ray} argument of
    {multifok_raytrace_make_frame}. It should trace a ray {R} towards an
    opaque scene or object and reports the first point {pHit(R)} of its surface
    that is hit by the ray, the normal direction {sNrm(R)} there, 
    and the color {sVal(R)}. It should return those results in
    of the hit point in {*pHit_P}, {*sNrm_P}, and {sVal[0..NC-1]}. 

    The ray will go through the point {pRay} with direction given by the
    unit-length vector {dRay}. Both {pRay}, {dRay}, and the returned hit
    point {pHit} are in some arbitary /scene coordinate system/.

    The ray is assumed to start at infinte distance from {pRay}
    in the direction opposite to {dRay}. The hit point may be anywhere
    along that ray, before or after {pRay}.
    
    If the procedure must always return a finite, non-{NAN} point in {*pHit_P}
    and some color in {sVal[0..NC-1]}.  It may return{(0,0,0)} in {*sNrm_P}
    to signify that the normal is undefined there.

    The {debug} parameter asks the procedure to print information
    about the ray tracing process. */

typedef void multifok_raytrace_img_to_scene_map_t(r2_t *p2_img, r3_t *p3_scene);
  /* Type of a procedure suitable for the {map_point} argument of
    {multifok_raytrace_make_frame}. It should map the 2D image coordinates
    {p2_img} of a point on the image plane to the corresponding
    3D scene coordinates {p3_scene}. */

typedef bool_t multifok_raytrace_debug_pred_t(i2_t *iPix);
  /* Type of a procedure suitable for the {debug_pixel} argument of
    {multifok_raytrace_make_frame}. It is given the row and colum
    indices of a pixel, and should return {TRUE} iff that
    procedure should print detailed information about that pixel. */

typedef void multifok_raytrace_report_ray_proc_t
  ( i2_t *iPix,
    i2_t *iSmp,
    double step,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay,
    double wRay,
    r3_t *pHit,
    double hHit,
    double vBlr,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  );
  /* Type of a procedure suitable for the {report_ray} argument of
    {multifok_raytrace_make_frame}. It is called for every ray cast for 
    the selected pixels.

    This procedure is given 
    
      * {iPix} row and column indices of the  pixel
      
      * {iSmp} indices of the sampoint relative to the pixel's center
      
      * {step} sampoint spacing (in pixels)
      
      * {wSmp} weight of that sampoint's contribution to the pixel
      
      * {pRay} sampoint's 3D scene system coordinates
      
      * {dRay} unit direction vector of the ray
      
      * {wRay} weight of that ray
      
      * {pHit} point where the ray hit the scene or object
      
      * {hHit} nominal height of the surface at that point
      
      * {vBlr} squared deviation of the ray from the ref direction {dRef}
      
      * {sNrm} surface normal at the hit point {pHit}
      
      * {sVal[0..NC-1]} apparent color of surface at that point
    
    The 2D coordinates of the sampoint will be {iPix.c[j] + 0.5 + step*iSmp.c[j]}
    for {j} in {0,1}.

    The parameters {pRay,dRay,pHit} will be in the scene's coordinate
    system. The image coordinates of the sub-sampling point {pSmp} may
    lie a bit outside the image domain {[0 _ NX] × [0 _ NY]}. */

typedef void multifok_raytrace_report_pixel_proc_t
  ( i2_t *iPix,
    r3_t *pCtr,
    double zFoc,
    double zDep,
    double shrp,
    double hAvg,
    double hDev,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  ); 
  /* Type of a procedure suitable for the {report_pixel} argument of
    {multifok_raytrace_make_frame}. It is called once for every pixel
    selected by {debug_pixel}.  This procedure is given:
    
      * {iPix} row and column indices of the pixel
      
      * {pCtr} center of the pixel in the 3D scene's coord sys
      
      * {zFoc,zDep} current focusing parameters
      
      * {shrp} average inferred sharpness parameter
      
      * {hAvg,hDev} average and deviation of scene height in pixel
      
      * {sNrm} average surface normal in pixel
      
      * {sVal[0..NC-1]} average pixel color */

void multifok_raytrace_paint_frame_rectangle
  ( multifok_frame_t *fr,
    int32_t xLo,
    int32_t xHi,
    int32_t yLo,
    int32_t yHi,
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    r3_t *dRef,
    double zFoc,
    double zDep,
    multifok_sampling_t *samp,
    bool_t verbose,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t *report_ray,
    multifok_raytrace_report_pixel_proc_t *report_pixel
  );
  /* Paints into a specific rectangle of the frame {fr} the result of 
    ray-tracing some generic object or scene with depth-of-focus blur.
    
    The frame's images {fr.sVal}, {fr.shrp}, {fr.hAvg}, {fr.hDev}, and
    {fr.sNrm} must have the same size. The painted rectangle consists of
    pixels in columns {xLo..xHi} and rows {yLo..yHi}. scene or object
    with depth-of-focus blur
    
    The value of each pixel {pix} in those images is computed by taking
    a set of /sampoints/ (pixel sub-sampling points) on the image plane
    around the pixel's center, converting each sampoint {pImg} to a
    point {pSmp} in scene coordinates with {map_point}, and tracing a
    number of rays through {pSmp}, as determined by the {samp} parameter
    record.  
    
    The {map_point} procedure should map points in the rectangle {[xLo _
    xHi+1] × [yLo _ yHi+1]} to the scene coordinates at {Z=zFoc} of the
    part of the scene covered by pixels {{xLo..xHi}×{yLo..yHi}} of the
    frame.
    
    The results of tracing these rays are combined as described under
    {multifok_raytrace_compute_point_properties} (q.v.) to produce the
    sampoint attributes {sVal(pSmp)}, {vBlr(pSmp)} and {hAvg(pSmp)}. The
    procedure combines these values of all {NS} sample points {pSmp} of
    the pixel with weights {wSmp[0..NS-1]} to obtain the values of
    {sVal(pix)}, {vBlr(pix)}, {hAvg(pix)}, {hDev(pix)}, and {sNrm(pix)} at that
    pixel. The value of {vBlr(pix)} is augmented to account for
    the spread of the sampoints relative to the pixel's center.

    The values of {sVal(pix)} are stored as the pixel values of the
    image {fr.sVal}. Ditto for the values of {hAvg(pix)}, {hDev(pix)},
    and {sNrm(pix)}. The value {vBlr(pix)} is converted to a sharpness
    indicator {shrp(pix) = 1/vBlr(pix)}, which is stored in the image
    {fr.shrp}.

    The fields {fr.zFoc} and {fr.zDep} are set to {zFoc} and {zDep}. The
    latter should be {+INF} for a sharp image (in which case {KR} should
    be 1). Otherwise the RMS spread of the ray directions will be one
    scene unit for each {zDep} of vertical distance from the sampling
    point.

    If {debug_pixel} is not null, it is called for every pixel in the
    image with the pixel indices {iPix}. If it returns true, the
    procedure prints detailed debugging information for the pixel with
    those indices. In that case, if {report_ray} is not {NULL}, the
    procedure also calls {report_ray} for every traced ray that
    contributed to that pixel's values. Also in that case, if
    {report_pixel} is not {NULL}, calls it once with the computed pixel
    properties. */

multifok_frame_t *multifok_raytrace_make_frame
  ( int32_t NC,
    int32_t NX,
    int32_t NY,
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    r3_t *dRef,
    double zFoc,
    double zDep,
    multifok_sampling_t *samp,
    bool_t verbose,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t *report_ray,
    multifok_raytrace_report_pixel_proc_t *report_pixel
  );
  /* Creates a frame record {fr} with images {fr.sVal}, {fr.shrp},
    {fr.hAvg}, {fr.hDev}, and {fr.sNrm}. All five images will have
    {NX} columns of {NY} pixels. The {sVal} image will have {NC}
    channels.  Then fills them with {multifok_raytrace_paint_frame_rectangle}
    with rectangle {{0..NX-1}×{0..NY-1}}. 
    
    The {map_point} procedure should map points in the rectangle {[0 _
    NX] × [0 _ NY]} to the scene coordinates at {Z=zFoc} of the
    part of the scene covered by the frame images. */

void multifok_raytrace_compute_point_properties
  ( multifok_raytrace_proc_t *trace_ray,
    r3_t *pObs, 
    uint32_t KR,
    r3_t *dRef,
    double zFoc,
    double zDep,
    bool_t debug,
    i2_t *iPix,
    i2_t *iSmp,
    double step,
    double wSmp,
    multifok_raytrace_report_ray_proc_t report_ray,
    r3_t *sNrm_pt_P,
    int32_t NC,
    float sVal_pt[],
    double *vBlr_pt_P,
    double *hAvg_pt_P,
    double *hVar_pt_P
  );
  /* Computes the mean normal {sNrm(pObs)}, color {sVal(pObs)}, blurring
    indicator {vBlr(pObs)}, height average {hAvg(pObs)}, and height
    variance {hVar(pObs)} seen from the point with scene coordinates
    {pObs}, by ray-tracing the scene or object with a set of {KR} rays
    through {pObs} and averaging the results. Returns these
    values in {*sNrm_pt_P}, {sVal_pt[0..NC-1]}, {*vBlr_pt_P},
    {*hAvg_pt_P}, and {*hVar_pt_P}.
    
    The direction {dir(R)} of each ray {R} will be chosen
    internally with a Gaussian distribution of deviations from the mean
    direction {dRef}. As explained under {multifok_raytrace_proc_t}, the
    rays will start an infinite distance away from {pObs} in the
    direction opposite to {dir(R)}; therefore, hits that occur
    ``behind'' {pObs} will be considered as well as those that occur
    ``ahead'' of it.
    
    As a special case, if {zDep} is {+oo}, then {KR} must be 1, and the
    procedure will shoot a single ray {R0} from the point {pObs}, with
    direction {dir(R0)=dRef}. This will result in "sharp" image, with
    no focus blur.
    
    The ray tracing is performed with
    {trace_ray(pRay,dRay,debug,&pHit,&sNrm,NC,sVal)} where {pRay=pObs} is the
    starting point, and {dRay=dir(R)} is the unit-length ray direction
    vector (both in scene coordinates). The result is an {NC}-channel
    color value {sVal(R)[0..NC-1]} and the scene coordinates of the hit
    point {pHit(R)}.
    
    The average normal {sNrm_pt} is computed by averaging the surface
    normals {sNrm(R)} at the hit points of all rays {R} and normalizing the
    result to unit length. If the result is invalid, it is replaced by
    {(0,0,0)}.

    The computed point color {sVal(pObs)}, its blurring indicator
    {vBlr(pObs)}, and average height {hAvg(pObs)} will be be the average
    of the ray colors {sVal(R}}, blurring indicator {vBlr(R)}, and
    {hHit(R)}. The variance {hVar(pObs)} is the variance of the values
    {hHit(R)} over all {R}. The value of {hDev(p)} is the deviation of
    {hVal(R)} over those rays.

    Let {vHit(R)} be the vector {pHit(R)-pObs}. The height {hHit(R)} will
    be {zFoc} minus the projection of {vHit(R)} along the direction {dRef},
    assumed to point perpendicular to the image plane and away from the
    camera. The tracing also produces a blurring indicator {vBlr(R)}
    which is the square of length of the projection of {vHit(R)}
    perpendicular to {dRef}.
    
    The burring indicator {vBlr(pObs)} of the point will be zero if
    {zDep} is {+INF} or the sampoint {pObs} lies on the scene's surface.
    In this latter case, {hAvg} will be {zFoc} and {hVar} will be zero.
    
    If {debug} is true, also prints debugging infrormation, and, if
    {report_ray} is not {NULL}, also calls it with the ray data.
    The parameters {iPix,iSmp,step,wSmp} are used only for this purpose. */

/* RAY-LEVEL LOGGING AND DEBUGGING */

void multifok_raytrace_write_ray_data
  ( FILE *wr,
    i2_t *iPix,
    double pixSize,
    i2_t *iSmp,
    double step,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay,
    double wRay,
    r3_t *pHit,
    double hHit,
    double vBlr,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  );
  /* Writes the data of a ray to file {wr}, in a format suitable for
    plotting and similar analyses. The parameters (other than {wr} and
    {pixSize}) are those of {report_ray} (q.v.). The line format is
    described in {multifok_raytrace_write_ray_data_INFO}. */

#define multifok_raytrace_write_ray_data_INFO \
  "The fields on each line are separated by blank space, with no" \
  " parentheses or other  delimiters.  They are:\n" \
  "\n" \
  "       {iPix.x} {iPix.y} the column and row indices of the pixel.\n" \
  "\n" \
  "       {pixSize} width and height of a pixel in scene coordinates.\n" \
  "\n" \
  "       {iSmp.x} {iSmp.y} the indices of the sampoint (sub-sampling point)" \
  " relative to the pixel's center, in some range {-HS..+HS}.\n" \
  "\n" \
  "       {step} the spacing between sampoints, in pixels.\n" \
  "\n" \
  "       {wSmp} the relative weight of that sub-sampling point" \
  " among those in the pixel.\n" \
  "\n" \
  "       {pRay.x} {pRay.y} {pRay.z} the scene coordinates {pRay} of that" \
  " sampling point.\n" \
  "\n" \
  "       {dRay.x} {dRay.y} {dRay.z} the scene coordinates of the ray's" \
  " unit direction vector {dRay}.\n" \
  "\n" \
  "       {wRay} the relative weight of that ray among those cast" \
  " through that sub-sampling point.\n" \
  "\n" \
  "       {pHit.x} {pHit.y} {pHit.z} the scene coordinates of the first" \
  " point {pHit} where the ray hit the scene's surface.\n" \
  "\n" \
  "       {hHit} the nominal height of the surface at that point.\n" \
  "\n" \
  "       {vBlr} the blurring indicator of that ray (square of" \
  " distance from {pRay} to {pHit} measured parallel to the image plane).\n" \
  "\n" \
  "       {sNrm.x} {sNrm.y} {sNrm.z} the surface normal at the point {pHit}.\n" \
  "\n" \
  "       {sVal[0]} {sVal[1]} {sVal[2]} the apparent surface color at the point {pHit}.\n" \
  "\n" \
  "    The unused color channels {sVal[NC..2]} are set to zeros."

void multifok_raytrace_show_ray_data
  ( FILE *wr,
    int32_t indent,
    i2_t *iPix,
    double pixSize,
    i2_t *iSmp, 
    double step,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay,
    double wRay,
    r3_t *pHit,
    double hHit,
    double vBlr,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  );
  /* Prints the data of a ray to file {wr}, in a legible format. The
    parameters (other than {wr} and {pixSize}) are those of {report_ray}
    (q.v.). */

/* PIXEL-LEVEL LOGGING AND DEBUGGING */

void multifok_raytrace_write_pixel_data
  ( FILE *wr,
    i2_t *iPix,
    double pixSize,
    r3_t *pCtr,
    double zFoc,
    double zDep,
    double shrp,
    double hAvg,
    double hDev,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  );
  /* Writes the computed data for pixel {oPix} to file {wr}, in a 
    format suitable for plotting
    and similar analyses.  The parameters (other than {wr})
    are those of {report_pixel} (q.v.).  The line format is
    described in {multifok_raytrace_write_pixel_data_INFO}. */

#define multifok_raytrace_write_pixel_data_INFO \
  "The fields on each line are separated by blank space, with no" \
  " parentheses or other  delimiters.  They are:\n" \
  "\n" \
  "       {iPix.x} {iPix.y} the column and row indices of the pixel.\n" \
  "\n" \
  "       {pixSize} width and height of a pixel in scene coordinates.\n" \
  "\n" \
  "       {pCtr.x} {pCtr.y} {pCtr.z} the scene coordinates of the pixel's" \
  " center {pCtr}.\n" \
  "\n" \
  "       {zFoc} the scene {Z} corrdinate of the simulated in-focus plane.\n" \
  "\n" \
  "       {zDep} the nominal depth of focus.\n" \
  "\n" \
  "       {shrp} the mena sharpness indicator, which is the reciprocal of square of" \
  " distance from the pixel sub-sampling point to the ray's" \
  " hit point, measured parallel to the image plane.\n" \
  "\n" \
  "       {hAvg} {hDev} the average and deviation of the {Z} coordinate" \
  " of the scene's hit points.\n" \
  "\n" \
  "       {sNrm.x} {sNrm.y} {sNrm.z} the average surface normal in the pixel.\n" \
  "\n" \
  "       {NC} The number of color channels used, in {1..3}.\n" \
  "\n" \
  "       {sVal[0]} {sVal[1]} {sVal[2]} The computed color of the pixel.\n" \
  "\n" \
  "    The averages {shrp,sNrm,sVal[0..NC-1],hAvg} and the deviation {hDev} are" \
  " computed over all the rays cast for the pixel, using the appropriate" \
  " weights.  The unused color channels {sVal[NC..2]} are set to zeros."

void multifok_raytrace_show_pixel_data
  ( FILE *wr,
    int32_t indent,
    i2_t *iPix,
    double pixSize,
    r3_t *pCtr,
    double zFoc,
    double zDep,
    double shrp,
    double hAvg,
    double hDev,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  );
  /* Prints the data of a pixel to file {wr}, in a legible format. The
    parameters (other than {wr} and {pixSize}) are those of {report_pixel}
    (q.v.). */

#endif
