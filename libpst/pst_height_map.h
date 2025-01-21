#ifndef pst_height_map_H
#define pst_height_map_H

/* pst_height_map.h -- procedures for working with height maps. */
/* Last edited on 2025-01-19 15:29:15 by stolfi */

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

float_image_t *pst_height_map_shrink(float_image_t *IZ);
  /* Given a height map {IZ}, with samples of a height function, returns
    another height map {JZ}, with half the size as {IZ} and the same
    number of channels, containing the samples of the height function at
    half the original resolution and half the original scale.
    
    The reduction is meant to be compatible with that of
    {pst_slope_map_shrink}.  Thus, if {IZ} has {NXI} cols and {NYI} rows, 
    the output map {JZ} must have {NXI/2+1} cols and {NY/2+1} rows,
    rounded DOWN.
    
    Moreover, a height sample in column {x} and row {y} of {JZ} should
    be the half of the weighted average of samples {IZ[x',y']} where
    {x'=2*x-dx}, {y'=2*y-dy}, and {dx,dy} range in {0..1}, if such
    pixels exist. The factor of 1/2 is meant to maintain the
    compatibility with {pst_slope_map_shrink}, that does not change the
    slope values even though the pixels are shrunk to half the size.
    
    If the map {IZ} has two channels, the weights for the averaging are
    taken from channel 1 of {IZ}. IN that case, the weight of the pixel
    in {JZ} is the minimum of the weights of the input samples used in
    the avergaing. If {IZ} has only one channel, the weigths are assumed
    to be all 1.0. */

float_image_t *pst_height_map_expand(float_image_t *JZ, int32_t NX, int32_t NY);
  /* The approximate inverse of {pst_height_map_shrink}.
    Given a height map {JZ}, magnifies it to produce a height
    map {IZ} twice as big.  The heights too are scaled by 2.
    
    The output map {IZ} will have {NX} columns and {NY} rows. The
    original image {JZ} must have size {NXI/2 + 1} by {NYI/2 + 1},
    rounded down. */

void pst_height_map_shift_to_zero_mean(float_image_t *Z);
  /* Subtracts the mean of all height values in {Z} from each height
    value, so that the mean height becomes zero. 
    
    The map {Z} must have 1 or 2 channels. If it has 2 channels, the
    mean is a weighted averagem, with the weights taken from channel 1.
    If {Z} has only one channel, the weights are assumed to be 1.0.
    However, if a height value is {NAN} or infinity, the corresponding
    weight is set to zero. */
  

/* COMPARISON */

float_image_t *pst_height_map_compare
  ( float_image_t *AZ,
    float_image_t *BZ,
    bool_t zero_mean,
    double *sAZP,
    double *sBZP,
    double *sEZP,
    double *sreP
  );
  /* Returns a height map {EZ} whose channel 0 is the difference {AZ-BZ} of the 
    two given height maps (which must have the same size). 
    If {zero_mean} is TRUE, subtracts the mean values so that
    the mean difference is zero.  Also sets the weight channel 1 of {EZ}
    to the product of the weight channels of {AZ} and {BZ}.
    
    Also returns in {*sAZP,*sBZP} the standard deviations {sZ,sBZ} of the values
    of {AZ,BZ}, each from its own mean value; in {*sEZP} the root-mean-square
    value {sEZ} of the error {EZ}; and in {*sreP} the relative error {sEZ/sMZ}
    where {sMZ = hypot(sAZ,sBZ)/sqrt(2)}. */
      
/* REPORTGING */

void pst_height_map_level_analyze_and_write
  ( char *filePrefix,
    int32_t level,
    int32_t iter,
    double change,
    float_image_t *CZ,
    float_image_t *RZ, 
    bool_t writeImages,
    bool_t writeError
  );
  /* A procedure for monitoring the progress of iterative
    integration, especially multiscale, and writing the solution.
    
    The image {CZ} must be a height map with 1 or 2 channels. The image
    {RZ}, if not NULL, is a reference height map, which must have 1 oe 2
    channels and the same size as {CZ}.
    
    If {RZ} is non-null, the procedure computes an image {EZ=CZ'-RZ'}
    where {CZ'} is {CZ} minus the average value of {CZ}; and ditto for
    {RZ'} from {RZ}. In the averages, each pixel is weighted by the
    product of the weights in the two maps.
    
    If {RZ} is non-null, the procedure also computes values
    {sZ=RMS(CZ')}, {sRZ=RMS(RZ')}, {sMZ=hypot(sZ,sRZ)/sqrt(2)},
    {sEZ=RMS(EZ)}, and {sre=sEZ/sMZ}.   In the the {RMS} values, the weight each
    squared sample error is the product of the two weight channels.
     
    If {writeImages} is true, the image {CZ} is written to disk file
    "{filePrefix}-{level}-{iter}-Z.fni". In that case, if {RZ} is not
    null, the procedure also writes the error image {EZ} to
    "{filePrefix}-{level}-{iter}-eZ.fni".

    If {writeError} is true and {RZ} is not null, the procedure also
    writes a file called "{filePrefix}-{level}-{iter}-eZ.txt"
    containing a single line single line
    
      "{level} {NX} {NY} {iter} {change} {sAZ} {sBZ} {sEZ} {sre}"
      
    If {level} is negative, the "-{level}" part is omitted from the file
    names. If {iter} is negative, the "-{iter}" part is omitted. If
    {level} is non-negative, all messages to {stderr} are indented by
    {2*level+2} spaces. */   

#endif
