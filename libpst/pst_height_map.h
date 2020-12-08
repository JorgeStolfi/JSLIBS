#ifndef pst_height_map_H
#define pst_height_map_H

/* pst_height_map.h -- procedures for working with height maps. */
/* Last edited on 2017-01-04 23:32:23 by stolfilocal */

#include <float_image.h>

/* HEIGHT MAPS
  
  A /height map/ is a single-channel float-valued image where
  the pixel value represents a height field {Z(X,Y)}. */

void pst_height_map_expand
  ( float_image_t *JZ,
    float_image_t *JW,
    float_image_t *IZ,
    float_image_t *IW
  );
  /* Given a height map {JZ} and its (possibly null) weight map {JW},
    magnifies it to produce a height map {IZ} twice as big,
    and (if {IW} is not null) a corresponding weight map {IW}.
    The heights too are scaled by 2.
    
    If the output map {IZ} has {NXI} columns and {NYI} rows, the
    original image {JZ} must have size {NXI/2 + 1} by {NYI/2 + 1},
    rounded down. If {JW} is given, it must have one channel and the
    same size as {JZ}, Ditto for {IW} and {IZ}. */

float_image_t *pst_height_map_shrink(float_image_t *IZ, int avgWidth);
  /* Given a height map {IZ}, with samples of a height function {Z},
    returns another height map {JZ}, with half the size as {IZ},
    containing the samples of the height function at half the original
    resolution and half the original scale. If the given image has
    size {NX} by {NY}, the result has size {NX/2+1} by {NY/2+1},
    rounded down. In particular, if {IZ} has {2^k+1} columns then {JZ}
    will have {2^{k-1}+1} columns.
    
    A sample in column {x} and row {y} of {JZ} is conceptually
    located at the corner {(x,y)} of {JZ}'s domain and at the corner {(2x,2y)} of {IZ}'s domain.
    The {avgWidth} parameter determines the width of the averaging kernel.
    If {avgWidth = 1}, uses subsampling (copies 1 every 2x2 pixels).
    If {avgWidth = 2}, averages a window of {3x3} pixels with binomial weights.
    If {avgWidth >= 3} , uses a smoother window of that width.
    
    Pixels of {IZ} that fall outside the domain are ignored in
    the averaging. */

/* 
  COMPARISON */

float_image_t *pst_height_map_compare
  ( float_image_t *AZ,
    float_image_t *BZ,
    float_image_t *U,
    bool_t zero_mean,
    double *sAZP,
    double *sBZP,
    double *sEZP,
    double *sreP
  );
  /* Returns a height map {EZ} which is the difference {AZ-BZ} of the 
    two given height maps (which msut have the same size). 
    If {zero_mean} is TRUE, subtracts the mean values so that
    the mean difference is zero. 
    
    Also returns in {*sAZP,*sBZP} the standard deviations {sZ,sBZ} of the values
    of {AZ,BZ}, each from its own mean value; in {*sEZP} the root-mean-square
    value {sEZ} of the error {EZ}; and in {*sreP} the relative error {sEZ/sMZ}
    where {sMZ = hypot(sAZ,sBZ)/sqrt(2)}. */
      
/* REPORTGING */

typedef void pst_height_map_report_proc_t(int level, int iter, double change, bool_t final, float_image_t *OZ); 
  /* Type of a client-given procedure that may be called
    by recursive integrators to report the height map used in each scale.
    The argument {iter} should be the number of iterations already done
    (0 = initial guess) and {change} should be the max height change from the 
    previous iteration (irrelevant for the initial guess).  The {final} arg should
    be true if the iteration has stopped. */   
   
void pst_height_map_level_analyze_and_write
  ( char *filePrefix,
    int level,
    bool_t levelTag,
    int iter,
    bool_t iterTag,
    double change,
    float_image_t *CZ,
    float_image_t *RZ,
    float_image_t *U, 
    bool_t writeImages,
    bool_t writeError,
    int indent
  );
  /* A procedure for monitoring the progress of iterative
    integration, especially multiscale, and writing the solution.
    
    The image {CZ} must be single-channel.
    
    The image {RZ}, if not NULL, is a reference height map, which
    must be single-channel and must have the same size as {CZ}.

    The image {U} should be a single-channel height confidence mask,
    with the same size as {CZ}. If null, it is assumed to be all 1's.
    
    If {RZ} is non-null, the procedure computes an image {EZ = CZ'
    - RZ'} where {CZ'} is {CZ} minus the average value of {CZ}; and
    ditto for {RZ'} from {RZ}. In the averages, each pixel is weighted
    by the corresponding pixel of {U}.
    
    If {RZ} is non-null, the procedure also computes values
    {sZ=RMS(CZ')}, {sRZ=RMS(RZ')}, {sMZ=hypot(sZ,sRZ)/sqrt(2)},
    {sEZ=RMS(EZ)}, and {sre=sEZ/sMZ}. All the {RMS} values use {U} as
    the weight mask.
     
    If {writeImages} is true, the images {CZ} is written to disk file
    "{filePrefix}-{level}-{iter}-Z.fni". In that case, if {RZ} is not
    null, the procedure also writes the error image {EZ} to
    "{filePrefix}-{level}-{iter}-eZ.fni".

    If {writeError} is true and {RZ} is not null, the procedure also
    writes a file called "{filePrefix}-{level}-{iter}-eZ.txt"
    containing a single line single line
    
      "{level} {NX} {NY} {iter} {change} {sAZ} {sBZ} {sEZ} {sre}"
      
    If {levelTag} is false the "-{level}" part is omitted from the file 
    names.  If {iterTag} is false the "-{iter}" tag is omitted.
    
    All messages to {stderr} are indented by {indent} spaces. */    

#endif
