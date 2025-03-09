#ifndef pst_height_map_H
#define pst_height_map_H

/* pst_height_map.h -- procedures for working with height maps. */
/* Last edited on 2025-03-02 12:10:15 by stolfi */

#include <float_image.h>

/* HEIGHT MAPS
  
  A /height map/ is a float-valued image with 1 or 2 channels, where
  the samples in channel 0 are height values {Z[x,y]}.  Channel 1,
  if present, stores the corresponding reliability weights {U[x,y]}.
  
  A height map {Z} is usually related to a slope (or normal) map {G}.
  In that situation, the PIXELS of {Z} correspond to CORNERS
  of pixels of {G}.  More precisely, pixel {Z[x,y]} is the low corner
  of pixel {G[x,y]}; that is it lies between pixels {G[x-dx,y-dy]} where
  {dx,xy} range in {0,1} -- when those pixels exist in {G}.  Conversely,
  {G[x,y]} is the pixel whose corners are {Z[x+dx,y+dy]} for {dx,dy} in {0,1}. */

/* MAP SHRINKING */

float_image_t *pst_height_map_shrink(float_image_t *IZ, double scale);
  /* Given a height map {IZ}, with samples of a height function, returns
    another height map {JZ}, with half the size of {IZ} and the same
    number of channels, containing the samples of the height function at
    half the original resolution. The height values (channel 0) of the
    result are multipled by the given {scale}.
    
    The reduction is meant to be compatible with that of
    {pst_slope_map_shrink}. Thus, if {IZ} has {NXI} cols and {NYI} rows,
    the output map {JZ} will have {NXJ=(NXI-1)/2+1} cols and {NYJ(NYI-1)/2+1}
    rows, rounded UP; that is, {NXJ=NXI/2+1} and {NYJ=NYI/2+1}, rounded
    DOWN. The number of channels of {JZ} will be the same as that of {IZ}
    (which must be 1 or 2).
    
    The {scale} factor is typically {0.5}, in order to maintain
    compatibility with {pst_slope_map_shrink}, that does not change the
    slope values even though the pixels are shrunk to half the size.
    
    The reduction is performed by {pst_vertex_map_shrink(IZ,1,NXJ,NYJ,scale)} 
    (q.v.). */

float_image_t *pst_height_map_expand(float_image_t *JZ, int32_t NXI, int32_t NYI, double scale);
  /* The approximate inverse of {pst_height_map_shrink}.
    Given a height map {JZ}, magnifies it to produce a height
    map {IZ} twice as big.  The height values (channel 0) of
    the result are multiplied by {scale}.
    
    The output map {IZ} will have {NXI} columns and {NYI} rows. The
    original image {JZ} must have {NXI/2+1} cols and {NYI/2+1} rows, rounded
    DOWN. The number of channels will be the same as that of {JZ} 
    (which must be 1 or 2).
    
    The {scale} factor is typically {2.0}, the inverse of 
    what should be used for {pst_height_map_shrink}.
    
    The expansion is performed by {pst_vertex_map_expand(JZ,1,NXI,NYI,scale)}
    (q.v.) */

/* PERTURBING */

void pst_height_map_perturb(float_image_t *A, int32_t wch, double relNoise, double absNoise);
  /* Adds to every sample of every channel {c} of {A}, except channel
    {wch}, a low-frequency perturbation in the range {[-mag _ +mag]},
    where {mag} is {hypot(absNoise,relNoise*(vmax-vmin))} and {[vmin _
    vmax]} is the original range of samples in that channel, ignoring
    non-finite samples. Uses {pst_height_map_perturb_channel(A.c.mag)}
    (q. v.).
    
    Sample {A[wch,X,Y]}, if it exists, is assumed to be a reliability
    weight (which must be finite and non-negative) for the other samples
    {A[c,X,Y]} of the same pixel. If any sample {A[c,X,Y]} was not
    finite, or the reliability weight of pixel {[X,Y]} is zero, then
    {A[wch,X,Y]} (if it exists) is set to zero, and all other samples
    of that pixels are set to {NAN}. */ 

void pst_height_map_perturb_channel(float_image_t *A, int32_t c, double mag);
  /* Adds to every sample of channel {c} of {A} a low-frequency
    perturbation in the range {[-mag _ +mag]}. Does not change samples
    that are not finite, or that overflow when the perturbation is
    applied. */
    
/* COMPARISON AND REPORTING */

void pst_height_map_analyze_and_write
  ( char *outPrefix,
    int32_t level,
    int32_t iter,
    double change,
    float_image_t *Z,
    float_image_t *R,
    bool_t verbose
  );
  /* A procedure for monitoring the progress of iterative
    integration, especially multiscale, and writing the solution.
    
    The image {Z} must be a height map with 1 or 2 channels. The image
    {R}, if not NULL, is a reference height map with same 
    channel count and the same size as {Z}.
    
    If {R} is non-null, the procedure calls {pst_height_map_compare} (q. v.)
    to computes an image {E=Z-R}, and the statistical parameters
    {devA,devB,undefE,avgE,devE,devRelE}. 
     
    The maps {Z}, and {E} are written to disk file
    "{iterPrefix}-{LL}-{NNNNNNNNN}-{tag}.fni" where {tag} is "Z", 
    or "E", respectively. Here {LL} is the {level} as "%02d" and
    {NNNNNNNNN} is {iter} as "%09d". However, if {iter} is 0 or negative
    the {NNNNNNNNN} part is replaced by "beg" and "end", respectively.
    The map {E} is written only if {R} is not {NULL}.
    
    If {R} is not {NULL}, the procedure also writes a file
    called "{iterPrefix}-{level}-{iter}-E.txt" containing a single line
    
      "{level} {NX} {NY} {iter} {change} {devA} {devB} {avgE} {devE} {devRelE}"

    If {level} is negative, the "-{level}" part is omitted from the file
    names. If {iter} is negative, the "{iter}" part is "fin"; if {iter} is zero,
    the {iter} in then. If
    {level} is non-negative, all messages to {stderr} are indented by
    {2*level+2} spaces. */   

#endif
