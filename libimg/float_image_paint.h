#ifndef float_image_paint_H
#define float_image_paint_H

/* Tools for drawing into float images. */
/* Last edited on 2023-11-26 06:42:18 by stolfi */ 

#include <bool.h>
#include <ellipse_crs.h>
#include <float_image.h>
#include <float_image_func.h>

/* IMAGE COORDINATES 

  For the functions in this interface, the domain of a float image {A}
  is an axis-aligned rectangle {[0 _ nx] × [0 _ ny]} of the plane
  {R^2} with pixel in row 0, column 0 adjacent to the origin, where
  {nx = A->sz[1]} and {ny = A->sz[2]}. Each pixel with indices {(x,y)}
  is conceptually a square with side 1, whose corners are {(x,y)} and
  {(x+1,y+1)}, and whose center is {(x+0.5,y+0.5)}.   */
  
/* RETURN VALUES
    
    The return value of all these procedures is the sum of the opacities
    of the overlaid samples (between 0 = transparent, 1 = opaque);
    namely, the value of the {mask} function, clipped to the range
    {[0_1]} and antialiased as above. This value is computed even if {A}
    is {NULL}. */
   
double float_image_paint_sample
  ( float_image_t *A, 
    int c, 
    int ix, 
    int iy, 
    float_image_func_t *func, 
    float_image_func_t *mask, 
    int m
  );
  /* Overlays pixel {ix,iy} of channel {c} of image {A} with the
    procedural image {func}, whose opacity is defined
    by the procedural image {mask}.
    
    The opacity {mask(x,y)} is interpreted as the local fraction of
    the background that is covered by the overlay. Thus,
    {mask(x,y)<=0} means that the overlay is transparent (preserving
    the underlying image), {mask(x,y)>=1} means it is fully opaque,
    and {mask(x,y)} between 0 and 1 means it is semi-transparent.
    
    As special cases, assumes that the overlay is fully transparent if
    {func(x,y) == NAN}, irrespective of {mask}; otherise, if {mask} is
    NULL, assumes {mask(x,y) == 1} for all {x,y}. Fails if {mask(x,y)}
    is {NAN}.
    
    The overlay and mask images are antialiased by averaging them over a grid of
    {2*m+1} by {2*m+1} subsamples within a window of size {2 × 2}
    pixels centered on pixel {ix,iy}, with C1 piecewise quadratic
    weights.  In particular, if {m} is zero, samples {func} and 
    {mask} only at the center of the pixel.  If {m} is 1, samples them
    at the pixel center, corners, and edge midpoints.  */

double float_image_paint_samples
  ( float_image_t *A, 
    int c, 
    int xLo, 
    int xHi,
    int yLo,
    int yHi, 
    float_image_func_t *func,
    float_image_func_t *mask,
    int m
  );
  /* Applies {float_image_paint_sample} to all pixels {ix,iy}
    with {ix} in {xLo..xHi} and {iy} in {yLo..yHi}. */

/* SYMBOL PAINTING

  The procedures in this section paint simple geometric 
  shapes on the given image {A}.  

  If the shape is closed, it is filled with the given sample value
  {vfill}, and then its outline is stroked with sample value {vdraw},
  using a round pen with radius {hwd} whose center is dragged around
  the outline.  
  
  An open shape can only be stroked, and admits no
  {vfill} parameter. 
  
  If {vfill} is {NAN}, the shape is not filled (so the original image
  values will remain). The outline is not drawn if {vdraw} is {NAN} or
  {hwd} is zero.

  The shape is antialiased by taking {2*m+1} subsamples per pixel
  along each axis, as per {float_image_paint_sample}.  */

double float_image_paint_rectangle
  ( float_image_t *A, 
    int c,           /* Channel. */                                
    double xmin,     /* Min X coordinate. */                  
    double xmax,     /* Max X coordinate. */                  
    double ymin,     /* Min Y coordinate. */                  
    double ymax,     /* Max Y coordinate. */                  
    double hwd,      /* Radius of pen tip. */  
    float vfill,     /* Ink value for filling. */                              
    float vdraw,     /* Ink value for stroking. */  
    int m            /* Subsampling parameter. */
  );
  /* Paints into channel {c} of image {A} the rectangle {[xmin_xmax]×[ymin_ymax]}.  
    The outline has sharp corners as if drawn with a square pen with 
    side {2*hwd}. The sign of {hwd} is ignored. */

double float_image_paint_dot
  ( float_image_t *A, 
    int c,           /* Channel. */                                
    double xctr,     /* Center's X coordinate. */                  
    double yctr,     /* Center's Y coordinate. */                  
    double rad,      /* Radius of cross (not counting {hwd}). */   
    double hwd,      /* Radius of pen tip. */                      
    bool_t round,    /* If {TRUE}, dor is round, else square. */
    bool_t diagonal, /* If {TRUE}, rotate 45 degrees. */           
    float vfill,     /* Ink value for filling. */                              
    float vdraw,     /* Ink value for stroking. */  
    int m            /* Subsampling parameter. */
  );
  /* Paints into channel {c} of image {A} a dot centered at {xctr,yctr}.
    
    If {round} is true, the dot is a circle with radius {rad} (pixels);
    otherwise it is it is a square with half-side {rad}. 
    In this second case, if {diagonal} is true, the square is rotated by 45 degrees.
    The signs of {rad} and {hwd} are ignored. */

double float_image_paint_smudge
  ( float_image_t *A, 
    int c,         /* Channel. */                                
    double xctr,   /* Center's X coordinate. */                  
    double yctr,   /* Center's Y coordinate. */                  
    double xdev,   /* Standard deviation in X direction. */   
    double ydev,   /* Standard deviation in Y direction. */                      
    float vfill,    /* Ink value at center. */  
    int m          /* Subsampling parameter. */
  );
  /* Paints onto channel {c} of image {A} a fuzzy dot centered at {xctr,yctr}.
    
    The dot will have color {vfill}, fully opaque at its
    center, fading out to transparency according to an axis-aligned
    Gaussian bell fuction with deviations {xdev,ydev} along each axis. */

double float_image_paint_cross
  ( float_image_t *A, 
    int c,           /* Channel. */
    double xctr,     /* Center's X coordinate. */
    double yctr,     /* Center's Y coordinate. */
    double rad,      /* Radius of cross (not counting {hwd}). */ 
    bool_t empty,    /* If {TRUE}, the center will be empty. */
    double hwd,      /* Radius of pen tip. */
    bool_t diagonal, /* If {TRUE}, rotate 45 degrees. */
    float vdraw,     /* Ink value for stroking. */  
    int m            /* Subsampling parameter. */
  );
  /* Draws into channel {c} of image {A} a cross centered
    at {xctr,yctr}, with arms of length {rad}. 
    
    If {empty} is true, the center of the cross will be left empty,
    and only the outer quarters of the arms will be drawn.  
    Note that the inner ends of the arms will have round caps that
    will partly fill the gap.
    
    If {diagonal} is true, the cross is rotated 45 degrees. 
    
    The signs of {rad} and {hwd} are ignored. */

/* ELLIPSE PAINTING */

double float_image_paint_ellipse_crs
  ( float_image_t *A,
    int c,            /* Channel. */
    ellipse_crs_t *E, /* Ellipse parameters. */
    double hwd,       /* Radius of pen tip. */
    float vfill,      /* Ink value for the interior, or {NAN}. */  
    float vdraw,      /* Ink value for stroking the outline, of {NAN}. */  
    int m             /* Subsampling parameter. */
  );
  /* Paints into channel {c} of {A} the image of ellipse {E}. */ 

double float_image_paint_ellipse_ouv
  ( float_image_t *A,
    int c,            /* Channel. */
    r2_t *ctr,        /* Center coords. */
    ellipse_ouv_t *F, /* Ellipse parameters, relative to center. */
    double hwd,       /* Radius of pen tip. */
    float vfill,      /* Ink value for the interior, or {NAN}. */  
    float vdraw,      /* Ink value for stroking the outline, of {NAN}. */  
    int m             /* Subsampling parameter. */
  );
  /* Paints into channel {c} of {A} the image of ellipse with center {ctr} and shape {F}. */ 

double float_image_paint_ellipse_aligned
  ( float_image_t *A,
    int c,            /* Channel. */
    r2_t *ctr,        /* Center coords. */
    r2_t *rad,        /* Radii in X and Y. */
    double hwd,       /* Radius of pen tip. */
    float vfill,      /* Ink value for the interior, or {NAN}. */  
    float vdraw,      /* Ink value for stroking the outline, of {NAN}. */  
    int m             /* Subsampling parameter. */
  );
  /* Paints into channel {c} of {A} the image of axis-aligned 
    ellipse with center {ctr} and semi-radii {rad.c[0],rad.c[1]}
    in X and Y, respectively. */ 

#endif
