#ifndef pst_signature_H
#define pst_signature_H

/* pst_signature.h -- procedures for computing normals from light signatures. */
/* Last edited on 2006-04-28 23:25:15 by stolfi */

#include <pst_basic.h>

#include <vec.h>

typedef struct signature_t
  { double *rin;   /* {sig[i] is relative intensity under light field {i}. */
    double mag;    /* Normalizing factor. */
    double var;    /* Total noise variance. */
  } signature_t;
  /* 
    A {signature_t} describes the behavior of a bit of surface under a
    set of {NF} of pseudo-monochromatic light fields, normalized to
    allow for fast collinearity matching.
    
    Specifically, {mag*rin[i]} is the apparent intensity of that bit
    of surface under light field {i}. The {rin} data is normalized
    per channel, that is, the sum of squares of {rin[0..N-1]} is 1.0.
    Thus {1/mag} is the normalizing factor used to convert raw
    intensities into normalizing intensities..
      
    The value of {var} is the total noise variance present in
    {rin[0..NF-1]}. */

signature_t pst_signature_new(int NF);
  /* Allocates a light signature for {NF} light fields. */

void pst_signature_print(FILE *wr, int NF, signature_t *sig);
  /* Prints on {wr} the light signature {sig}, assumed 
    to be for {NF} light fields. */

void pst_signature_extract
  ( image_vec_t *IMGV,    /* Scene images under {NF} different light fields. */
    int maxval,        /* Maxval of original (quantized) images. */                 
    double noise,      /* Additional per-sample noise in images. */                 
    int c,             /* Channel. */                      
    int x,             /* Pixel column. */                                           
    int y,             /* Pixel row (0 = bottom). */  
    signature_t *sig   /* (OUT) Light signature. */
  );
  /* Extracts from channel {c} of the images {IMGV[0..NF-1]} (where
    {NF=IMGV.ne}) the light signature of the pixel in column {x} and row
    {y}. The signature is returned in {*sig}, which must have been
    pre-allocated for {NF} light fields.
    
    The total noise variance {sig->var} is estimated by assuming that
    each image sample was originally quantized with the given {maxval}
    and had extra noise with standard deviation {noise} added to it.
    As a special case, if {maxval} is 0, the input samples are assumed
    to be free from quantization errors. */

vec_typedef(signature_vec_t,signature_vec,signature_t);
  /* A vector of normalized light signatures. */

typedef struct light_table_t /* A table that maps light signatures to normal vectors. */
  { r2_t pos;              /* Position in scene images where table is most valid. */
    int NF;                /* Number of light fields. */             
    int NC;                /* Number of color channels. */
    int NE;                /* Number of color channels. */
    r3_vec_t nrm;          /* {nrm[k]} is the normal vector associated with entry {k}. */
    signature_vec_t *sig;  /* {sig[c][k]} is the light signature in channel {c} for that normal. */
  } light_table_t;
  /* A light table stored the normalized light signatures extracted
    from the light gauge images. Each entry {k} of the table comes from
    some pixel of those images, that was considered valid by some
    criterion.  The table contains the normal at that pixel, and
    the normalized light signatures for each color channel. */

light_table_t *pst_signature_build_table
  ( r2_t *pos,            /* Position in scene images where table is most valid. */
    float_image_t *NRM,   /* Normal map of light gauge. */ 
    image_vec_t *IMGV,    /* Images of the gauge object under various lightings. */ 
    bool_t cubify         /* TRUE uses the cube acceleration. */
  );
  /* Creates a light table from images {IMGV[0..NF-1]} of a light gauge under
    {NF} different light fields, where {NF = IMGV.ne}. Assumes that
    {NRM} is the gauge's /normal map/, that is, a three-channel image
    where the value of each pixel is an outwards-pointing unit vector
    perpendicular to the gauge surface, averaged over that pixel.
    Pixels where {NRM} is (0,0,0) are ignored. 
    
    The images {IMGV[0..NF-1]} must have the same dimensions and the same 
    number {NC} of channels. If {cubify} is true, the table will be
    provided with a bucket grid accelerator. */

void pst_signature_search_table
  ( int NF,               /* Number of light fields. */
    int NC,               /* Number of color channels. */
    signature_t sig[],    /* Normalized light signature for each channel (size {NC}). */ 
    r2_t pos,             /* Nominal position of {sig} in scene images. */
    int NG,               /* Number of light-to-normal tables. */
    light_table_t *tab,   /* Light-to-normal table extracted from gauge. */
    double *dsq,          /* (OUT) Discrepancy squared between {sig} and best match from table. */
    r3_t *nrm,            /* (OUT) Normal associated to best match in table. */
    float clr[]           /* (OUT) Intrinisc scene color (size {NC}). */
  );
  /* Scans the table {tab} for the best match of light signatures
    {sig[0..NC-1]}. Returns in {*diff} the expected distance
    (Euclidean) between {sig} and the best match from the tables.
    Stores in {*nrm} the normal vector associated to the best match,
    and in {clr[c]} the inferred intrinsic lightness of the scene's
    surface in channel {c}. 
    
    The procedure also stores in {dsq} the squared distance between
    {sig} and the best-matching signature from the table. The value of
    {dsq=0} means a perfect match; the value {dsq=1} means a complete
    mismatch. */

void pst_signature_normals_from_photos
  ( int NG,                  /* Number of light gauges. */
    light_table_t *tab[],    /* {tab[0..NG-1]} are the light-to-normal tables of the gauges. */
    image_vec_t *IMGV,       /* {IMGV[0..NF-1]} are photos of a scene under {NF} light fields. */
    int maxval,              /* Number of quantization levels in original quantized images. */
    double noise,            /* Standard deviation of additional per-sample noise. */
    float_image_t *NRM,      /* (OUT) Nomal map of the scene. */
    float_image_t *CLR,      /* (OUT) Intrinsic color map of the scene. */
    float_image_t *DIF       /* (OUT) Strangeness map of the scene. */
  );
  /* Produces a normal map {NRM}, an intrinsic color map {CLR}, and
    a strangeness map {DIF} from a list {IMGV[0..NF-1]} of photos of a
    scene under {NF} different light fields (where {NF=IMGV.ne}).
    
    The normal {NRM[p]} ate ach pixel {p} is estimated by comparing
    the appearance of {p} in the photos {IMGV[0..NF-1]} with the
    appearances of pixels from the light gauges, stored in the
    light-to-normal tables {tab[0..NG-1]}. The procedure uses the
    table {tab[g]} whose nominal positon {tab[g]->pos} is closest to
    {p}'s center. The procedure scans that table for a gauge pixel {q}
    whose appearance, under those {NF} light fields, best matches
    those of the scene pixel, except for possible differences in
    intrinsic color. The surface orientation at the {p} is then
    assumed to be the same as that of {q}. In the normal image, the
    three channels are the X, Y, and Z component of the
    outwards-pointint unit-length vector normal to the surface.
    
    The intrinsic color {CLR[p]} of the scene at each pixel {p} is
    obtained by comparing the overall intensity of that pixel in each
    color channel {c} with the overall intensity of the gauge's
    best-matching pixel {q} in that same channel. In this image, a
    sample value of 1 means that the scene is as reflective as the
    gauge in that particular channel.
    
    The strangeness {DIF[p]} of each pixel {p} is a measure of the
    difference between its appeareance in the photos and the
    appearance of the matching gauge pixel, after accounting for
    differences in intrinsic color. The value 0 means a perfect match;
    the value 1 means the largest mismatch possible. */


#endif
