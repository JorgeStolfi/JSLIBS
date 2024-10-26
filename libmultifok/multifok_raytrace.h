/* Generic raytracing for the {multifok} library. */
/* Last edited on 2024-10-22 12:15:16 by stolfi */

#ifndef multifok_raytrace_H
#define multifok_raytrace_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <i2.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>

#include <multifok_frame.h>

#define ix_DEBUG 10
#define iy_DEBUG 17
  /* Default pixel indices for which debugging info is to be printed out. */

typedef void multifok_raytrace_proc_t
  ( r3_t *pRay,
    r3_t *dRay,
    bool_t debug,
    r3_t *pHit_P,
    int32_t NC,
    float colr[]
  );
  /* Type of a procedure suitable for the {trace_ray} argument of
    {multifok_raytrace_make_frame}. It should trace a ray {R} towards an
    opaque scene or object and reports the first point of its surface
    that is hit by the ray. It should return the coordinates {pHit(R)}
    of the hit point, and the local color coordinates {colr(R)} as seen
    from the ray's source. These results should be returned in
    {*pHit_P}, and {colr[0..NC-1]}, respectively.

    The ray will go through the point {pRay} with direction given by the
    unit-length vector {dRay}. Both {pRay}, {dRay}, and the returned hit
    point {pHit} are in some arbitary /scene coordinate system/.

    The ray is assumed to start at infinte distance from {pRay}
    in the direction opposite to {dRay}. The hit point may be anywhere
    along that ray, before or after {pRay}.

    If the ray fails to hit anything, it should return {(NAN,NAN,NAN)}
    in {*pHit_P} and all {NAN} in {colr[0..NC-1]}.

    The {debug} parameter asks the procedure to print information
    about the ray tracing process. */

typedef void multifok_raytrace_img_to_scene_map_t(r2_t *p2_img, r3_t *p3_scene);
  /* Type of a procedure suitable for the {map_point} argument of
    {multifok_raytrace_make_frame}. It should map the 2D image coordinates
    {p2_img} of a point on the image plane to the corresponding
    3D scene coordinates {p3_scene}. */

typedef bool_t multifok_raytrace_debug_pred_t(i2_t *iPix);
  /* Type of a procedure suitable for the {debug_pix} argument of
    {multifok_raytrace_make_frame}. It is given the row and colum
    indices of a pixel, and should return {TRUE} iff that
    procedure should print detailed information about that pixel. */

typedef void multifok_raytrace_report_ray_proc_t
  ( i2_t *iPix,
    r2_t *pSmp,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay,
    double wRay,
    r3_t *pHit,
    double hHit,
    double vBlr
  );
  /* Type of a procedure suitable for the {report_ray} argument of
    {multifok_raytrace_make_frame}. It is called for every ray that
    contributes to selected pixels.

    This procedure is given the row and column indices {iPix} of the
    pixel, the 2D image system coordinates {pSmp} of a subsampling point
    that contributed to that pixel, the weight {wSmp} of that point's
    contribution to the pixel, its 3D scene system coordinates {pRay},
    the unit direction vector {dRay} of a ray that contributed to that
    sub-sampling point's values, the weight {wRay} of that contribution,
    the point {pHit} where the ray through {pRay} with direction {dRay}
    hit the scene or object, the nominal height {hHit} of the surface at
    that point, and the squared deviation {vBlr} of the ray from the
    reference direction {dRef}.

    The parameters {pRay,dRay,pHit} will be in the scene's coordinate
    system. The image coordinates of the sub-sampling point {pSmp} may
    lie a bit outside the image domain {[0 _ NX] × [0 _ NY]}. */

multifok_frame_t *multifok_raytrace_make_frame
  ( int32_t NC,
    int32_t NX,
    int32_t NY,
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    r3_t *dRef,
    double zFoc,
    double zDep,
    int32_t NS,
    r2_t uSmp[],
    double wSmp[],
    int32_t NR,
    r2_t tRay[],
    double wRay[],
    bool_t verbose,
    multifok_raytrace_debug_pred_t *debug_pix,
    multifok_raytrace_report_ray_proc_t *report_ray
  );
  /* Creates a frame record {fr} with images {fr.sVal}, {fr.shrp}, {fr.hAvg},
    and {fr.hDev} by ray-tracing some generic scene or object with
    depth-of-focus blur. All four images will have {NX} columns of {NY} pixels.
    The {sVal} image will have {NC} channels.

    The value of each pixel {pix} of the four images is computed by
    taking a set of {NS} sample points on the image plane, displaced by
    {uSmp[0..NS-1]} from the center of {pix}. The procedure then traces
    a total of {NR} rays through these sample points with directions
    defined by the relative ray tilts {tRay[0..NR-1]} and the nominal
    depth-of-focus {zDep}. The number of rays {NR} must be an integer
    multiple of the number of sampling points {NS}; for each sampling
    point, the procedure will cast a disjoint subset of the rays with size 
    {KR=NR/NS}.  However, if {NR} is 1, it will use that single ray,
    assumed to be straight down ({tRay[0]=(0,0)}) for all sampling points,
    so {KR} will be 1.

    The result of tracing each ray {R} is an {NC}-channel color value
    {colr(R)[0..NC-1]} and the scene coordinates of the hit point
    {pHit(R)}. The ray tracing is performed with
    {trace_ray(pRay,dRay,debug,&pHit,NC,colr)} where {pRay =
    map_point(pSmp)} is the sampling point, {dRay} is the unit-length
    ray direction vector (both in scene coordinates), and {debug} is the
    result returned by {debug_pix(pix)}.

    Let {uRay} be the vector {pHit(R)-pSmp}. The height {zVal(R)} will
    be {zFoc} minus the projection of {uRay} along the direction {dRef},
    assumed to point perpendicular to the image plane and away from the
    camera. The tracing also produces a blurring indicator {vBlr(R)}
    which is the square of length of the projection of {uRay}
    perpendicular to {dRef}.

    The values {colr(pSmp)}, {vBlr(pSmp)} and {hAvg(pSmp)} are the
    averages of {colr(R)}, {vBlr(r)}, and {zVal(R)} over all the {KR}
    rays assigned to that sampling point, with the corresponding weights
    from {wRay[0..NR-1]}. The value of {Zdev(p)} is the weighted deviation
    of {zVal(R)} over those rays.

    Finally, the procedure combines these values of all {NS} sample points
    {pSmp} of the pixel with weights {wSmp[0..NS-1]} to obtain
    the values of {colr(pix)}, {vBlr(pix)}, {hAvg(pix)}, and {hDev(pix)}
    at that pixel. The value of {vBlr(pix)} is augmented with 1.0 to
    account for the fact that antialiasing (subsampling and averaging)
    causes some blurring.

    The values of {colr(pix)}, {hAvg(pix)}, and {hDev(pix)} are stored
    as the pixel values of the image {fr.sVal}. The value {vBlr(pix)} is
    converted to a sharpness indicator {shrp(pix) = 1/vBlr(pix)}, which
    is stored in the image {fr.shrp}.

    The fields {fr.zFoc} and {fr.zDep} are set to {zFoc} and {zDep}. The
    latter should be {+INF} if {NR} is 1, otherwise it should be inversely
    proportional to the mean spread of the ray directions {dRay}.

    For best results, both the ray tilts {tRay[0..NR-1]} and the
    sampling point displacements {uSmp[0..NS-1]} should be sorted
    by increasing distance from the origin {(0,0)}.

    If {debug_pixel} is not null, it is called for every pixel in the image
    with the pixel indices {iPix}. If it returns true, the procedure prints
    detailed debugging information for the pixel with those indices. In
    that case, if {report_ray} is not {NULL}, the procedure also calls
    {report_ray} for every traced ray that contributed to that pixel's
    values. */

r3_t multifok_raytrace_compute_ray_direction(r2_t *tRay, double aRay, double zDep);
  /* Returns a unit ray direction vectors {dRay} to
    ray-trace a scene from some image sampling point, given the
    relative ray tilt {tRay} and a ray pattern rotation angle {aray}.

    The ray will be pointing down ({dRay.c[2] < 0}).

    The ray will hit a horizontal plane located {zDep/2} below its
    starting point {p} at a point displaced from the vertical
    through {p} by the 2-vector {tRay} rotated by {aRay} radians
    about that vertical. */

void multifok_raytrace_show_ray_data
  ( FILE *wr,
    i2_t *iPix,
    double pixSize,
    r2_t *pSmp,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay,
    double wRay,
    r3_t *pHit,
    double hHit,
    double vBlr
  );
  /* Prints the data of a ray to file {wr}, in a legible format. The
    parameters (other than {wr}) are those of {report_ray} (q.v.). */

void multifok_raytrace_write_ray_data
  ( FILE *wr,
    i2_t *iPix,
    double pixSize,
    r2_t *pSmp,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay,
    double wRay,
    r3_t *pHit,
    double hHit,
    double vBlr
  );
  /* Writes the data of a ray to file {wr}, in a format suitable for plotting
    and similar analyses.  The parameters (other than {wr})
    are those of {report_ray} (q.v.).  The line format is
    described in {multifok_raytrace_write_ray_data_INFO}. */

#define multifok_raytrace_write_ray_data_INFO \
  "The fields on each line are separated by blank space, with no" \
  " parentheses or other  delimiters.  They are:\n" \
  "\n" \
  "       {iPix.x} {iPix.y} the column and row indices of the pixel.\n" \
  "\n" \
  "       {pixSize} width and height of a pixel in scene coordinates.\n" \
  "\n" \
  "       {pSmp.x} {pSmp.y} the image coordinates of the pixel" \
  " sub-sampling point {pSmp}.\n" \
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
  " distance from {pRay} to {pHit} measured parallel to the image plane)."

#endif
