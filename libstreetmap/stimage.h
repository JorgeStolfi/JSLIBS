/* stimage - images that represent geographic data. */
/* Last edited on 2003-08-26 19:58:56 by stolfi */

#ifndef stimage_H
#define stimage_H

#include <stmap.h>

typedef struct DataImage
  { int mmXLo;  /* Lowest X sample coordinate, in mm. */
    int nx;     /* Number of samples in X direction. */
    int mmYLo;  /* Lowest Y sample coordinate, in mm. */
    int ny;     /* Number of samples in Y direction. */
    int mmStep; /* Step between samples, in mm. */
    double *d;  /* Sample values. */
  } DataImage;
  /* A real-valued function of position, defined by interpolation
    on a rectangular array of sampled values {d[0..nx-1,0..ny-1]}. 
    Sample {[ix,iy]} is located at point 
    {(mmXLo + ix*mmStep, mmYLo + iy*mmStep)/1000} 
    and is stored in {d[ix + nx*iy]}. */

DataImage *st_image_new
  ( int mmXLo,  /* Lowest X sample coordinate, in mm. */
    int nx,     /* Number of samples in X direction. */
    int mmYLo,  /* Lowest Y sample coordinate, in mm. */
    int ny,     /* Number of samples in Y direction. */
    int mmStep  /* Step between samples, in mm */
  );
  /* Allocates a new {DataImage} with given attributes.
     The sample values are set to zero. */

DataImage *st_image_for_rect(Interval rx, Interval ry, double step);
  /* Allocates a new {DataImage} that spans the rectangle {rx × ry}
    (in meters) with sample points spaced by approximately {step}
    meters. The sample values are set to zero. */

int st_image_locate(DataImage *img, double x, double y);
  /* Returns the index into {img->d} of the sample closest 
    to point {(x,y)} (in meters), or -1 if that point lies 
    outside the rectangle covered by the samples. */
    
double st_image_get(DataImage *img, double x, double y);
  /* Returns the value of the sample of {img} closest to point {(x,y)}
    (in meters), or 0 if that point lies outside the
    rectangle covered by the samples. */

void st_image_set(DataImage *img, double x, double y, double val);
  /* Stores {val} into the sample of {img} closest to point {(x,y)}
    (in meters).  A no-op if that point lies outside the
    rectangle covered by the samples. */
    
void st_image_increment(DataImage *img, double x, double y, double val);
  /* Like {st_image_set}, but adds {val} to the sample instead of 
    storing it. */
    
void st_image_increment(DataImage *img, double x, double y, double val);
  /* Increments by {val} the pixel of {img} closest to {(x,y)}
    (in meters). */

double st_image_interpolate(DataImage *img, double x, double y);
  /* Returns the value of {img} at the point {(x,y)} (in meters),
    obtained by interpolating linearly between the four nearest
    samples. Samples outside the covered rectangle are assumed to be
    zero. */

Interval st_image_sample_range(DataImage *img);
  /* Returns the minimum and maximum sample values in image {img},
    as an interval. */

DataImage *st_image_blur(DataImage *img, double rad);
  /* Returns the convolution of the image {img} by a Gaussian 
    hump of total variance {rad^2}. */

void st_image_write_as_pgm
  ( FILE *wr, 
    DataImage *img,
    double bval,
    double wval
  );
  /* Writes the data image {img} to file {wr}, as a PGM image, mapping
    {bval} to black and {wval} to white. The attributes
    {mmXLo,mmYLo,mmStep,bval,wval} are also written at the end of the
    file, after the last pixel. */

DataImage *st_image_read_from_pgm(FILE *rd);
  /* Reads a data image from the given PGM file, which must have been
    written with {st_image_write_as_pgm}. In particular, expects to find the attributes
    {mmXLo,mmYLo,mmStep,bval,wval} at the end of the file. */

#endif
