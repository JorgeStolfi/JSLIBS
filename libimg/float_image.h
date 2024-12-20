#ifndef float_image_H
#define float_image_H

/* Multichannel images with floating-point samples. */
/* Last edited on 2024-12-20 17:49:36 by stolfi */ 

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <ix.h>

#define INF INFINITY

typedef struct float_image_t 
  { ix_size_t sz[3]; /* Image sizes: sz[0] = channels, sz[1] = columns, sz[2] = rows. */
    ix_step_t st[3]; /* Indexing step along each axis. */
    ix_pos_t bp;     /* Position of sample {[0,0,0]}. */
    float *sample;   /* Linearized sample vector. */
  } float_image_t; 
  /* The sample value for channel {c}, column {x}, and row {y} 
    of the image is stored in {sample[bp + c*st[0] + x*st[1] + y*st[2]]}.
    
    In the comments that follow, we use {NC}, {NX}, and {NY} to mean
    the number of channels, columns, and rows of the relevant image, 
    respectively. */

#define float_image_max_size 65356
  /* Maximum dimension along any axis (2^16), for sanity checks. */

/* ALLOCATION */

float_image_t *float_image_new (int32_t NC, int32_t NX, int32_t NY);
  /* Allocates a new multichannel image with {NC} channels,
    {NX} columns, and {NY} rows. */

float_image_t *float_image_copy (float_image_t *A);
  /* Creates a new float image with the same size and number
    of channels as {A}, and copies into its samples the samples of {A}. */

float_image_t *float_image_crop
  ( float_image_t *A, 
    int32_t cLo,
    int32_t cHi,
    int32_t xLo,
    int32_t xHi,
    int32_t yLo,
    int32_t yHi,
    float bg
  );
  /* Creates a new float image containing a copy of channels
    {cLo..cHi-1} of {A} clipped to the rectangle
    {[xLo_xHi]×[yLo_yHi]} (that is, to columns {xLo..xHi-1} and
    rows {yLo..yHi-1}). Substitutes {bg} for any samples that 
    do not exist in {A}. */

void float_image_free(float_image_t *A);
  /* Discards the pixel array *and the header* of the image {A}. */

/* ACESSING INDIVIDUAL SAMPLES */

float float_image_get_sample(float_image_t *A, int32_t c, int32_t x, int32_t y);
  /* Returns the sample in channel {c}, column {x}, row {y} of image {A}.
    Fails if any of {c,x,y} are out of bounds. */

float *float_image_get_sample_address(float_image_t *A, int32_t c, int32_t x, int32_t y);
  /* Returns the address of the sample in channel {c}, column {x}, row {y} 
    of image {A}. Fails if any of {c,x,y} are out of bounds. */

void float_image_set_sample(float_image_t *A, int32_t c, int32_t x, int32_t y, float v);
  /* Stores {v} into the sample in channel {c}, column {x}, row {y} of image {A}.
    Fails if any of {c,x,y} are out of bounds. */

/* ACESSING WHOLE PIXELS */

void float_image_get_pixel(float_image_t *A, int32_t x, int32_t y, float v[]);
  /* Copies into {v[c]} the sample in channel {c}, column {x} and row
    {y} of image {A}, for {c=0..NC-1}; where {NC} is the number of
    channels of {A}. Fails if any of {x,y} are out of bounds. */

void float_image_set_pixel(float_image_t *A, int32_t x, int32_t y, float v[]);
  /* Stores the sample value {v[c]}, into channel {c}, column {x} and
    row {y} of image {A}, for {c=0..NC-1}; where {NC} is the number
    of channels of {A}. Fails if any of {x,y} are out of bounds. */
    
void float_image_fill_pixel(float_image_t *A, int32_t x, int32_t y, float v);
  /* Stores the sample value {v} into channel {c}, column {x} and row
    {y} of image {A}, for {c=0..NC-1}; where {NC} is the number of
    channels of {A}. Fails if any of {x,y} are out of bounds. */

/* WINDOW EXTRACTION */

void float_image_get_window_samples
  ( float_image_t *A, 
    int32_t c, 
    int32_t x, 
    int32_t y, 
    int32_t nwx, 
    int32_t nwy,
    bool_t rep, 
    float v[]
  );
  /* Copies into {v[k]} the samples of channel {c} of {A} that lie
    inside a {nwx} by {nwy} window centered on column {x} and row {y} of
    image {A}. The widths {nwx,nwy} must be odd.
    
    The sample in column {ix} and row {iy} of the window is stored into
    {v[ix + nwy*iy]}.
    
    Fails if {c} is out of bounds. Any samples that fall outside the
    bounds of {A} are set either to {NAN} (if {rep} is false) or to the
    nearest in-bounds sample of channel {c} (if {rep} is true). */

void float_image_get_window_pixels
  ( float_image_t *A, 
    int32_t x, 
    int32_t y, 
    int32_t nwx, 
    int32_t nwy, 
    bool_t rep, 
    float v[] 
  );
  /* Copies into {v[k]} the pixels of {A} that lie inside a {nwx} by
    {nwy} window centered on column {x} and row {y} of image {A}. The
    widths {nwx,nwy} must be odd.
    
    The sample from channel {ic} of the image in column {ix} and row
    {iy} of the window is stored into {v[ic + NC*(ix + nwy*iy)]}.
    
    Any samples that fall outside the bounds of {A} are set either to
    {NAN} (if {rep} is false) or to the nearest in-bounds pixel (if
    {rep} is true). */

/* LOCAL MEAN AND VARIANCE */

void float_image_get_local_avg_var
  ( float_image_t *A, 
    int32_t c, 
    int32_t x, 
    int32_t y,
    int32_t hw,
    double wt[],
    double *avgP, 
    double *varP
  );
  /* Computes the local weighted average {*avgP} and local weighted
    variance {*varP} around pixel {x,y} of channel {c} of {A}.
    
    Uses a weight mask with {2*hw+1} columns and rows, centered over
    pixel {x,y} of {A}, with weight profile {wt[0..2*hw]}.
    
    Specifically, the weight of a sample in column {x+dx} and row
    {y+dy} ofgrep {A} is {wt[dx+hw]*wt[dy+hw]}, for {dx,dy} in {-hw..+hw};
    except that the weight is 0 if that sample is not in {A},
    or is {NAN}.
    
    The indices {x,y} need not be inside the domain of {A},
    but if all weights are zero (that is, if either {x} or {y} is 
    outside {A}'s domain by {hw} or more pixels) then the result is 
    NAN. */

void float_image_local_avg_var
  ( float_image_t *A, 
    int32_t cA, 
    int32_t hw,
    double wt[],
    float_image_t *M, 
    int32_t cM,
    float_image_t *V, 
    int32_t cV
  );
  /* Computes the local weighted average and local weighted variance
    around each pixel of channel {cA} of {A}, using
    {float_image_get_local_avg_var} with a weight mask of {2*hw+1}
    columns and rows and weight profile {wt[0..2*hw]}.
    
    If {M} is not null, stores the local averages into channel {cM} of
    {M}. If {V} is not null, stores the local variances into channel
    {cV} of {V}. These destination channels must be storage-disjoint
    from channel {cA} of {A}.
    
    If {M}'s width is different from that of {A}, the difference {Dx}
    of the widhs must be even, and then {M} is conceptually shifted by
    {Dx/2} pixels so that the centers of the domains are aligned.
    Ditto for the height of {M} and the dimensions of {V}. In any case
    the averages are computed over the min enclosing box of all three
    domains.
    
    This procedure may be used for smoothing. Although it is more
    expensive than FFT-based filtering, it does not assume that
    the image is periodic, and therefore does not suffer from
    wrap-around spillover. */

/* IMAGE GRADIENT */

void float_image_get_gradient_sobel
  ( float_image_t *A, 
    int32_t c, 
    int32_t x, 
    int32_t y,
    double *fxP, 
    double *fyP
  );
  /* Computes the gradient {*fxP,*fyP} at the center of pixel {x,y}
    of channel {c} of {A}. 
    
    Uses the Sobel operator that performs smoothing in the direction
    transversal to the differencing direction.  The X derivative {*fxP}
    will be 0 along the left and right edges, and the Y derivative 
    {*fyP} will be 0 along the top and bottom edges. */

/* IMAGE DILATION AND EROSION */

double float_image_get_dilated
  ( float_image_t *A, 
    int32_t c, 
    int32_t x, 
    int32_t y,
    int32_t hw,
    double wt[]
  );
  /* Returns the dilation of channel {c} of image {A} at pixel {x,y}. Uses a
    weight mask with {2*hw+1} columns and rows, with weights
    {wt[dx+hw]*wt[dy+hw]} for {dx,dy} in {-hw..+hw}. 
    
    More precisely, computes
    
      {max { wt[dx+hw]*wt[dy+hw]*A[c,x+dx,y+dy] : dx,dy \in -hw..+hw }}
    
    where {A[c,x,y]} is {float_image_get_sample(A,c,x,y)}. The vector
    {wt} must be a unidimensional weight table with {2*hw+1} entries.
    Usually the entries are non-negative and symmetric around
    {wt[hw]}. */

/* !!! Add {float_image_get_eroded} !!! */

/* ACESSING WHOLE ROWS OF SAMPLES */

void float_image_get_sample_row
  ( float_image_t *A, 
    int32_t c, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t y, 
    float v[] 
  );
  /* Sets {v[x-xmin]} to the sample of channel {c}, column {x}, row {y}
    of image {A}, for {x} in {xmin..xmax}; or to {NAN} if {x} is not in
    {0..NX-1} -- where {NX} is the number of columns of {A}.
    
    If {xmin>xmax} no samples are copied. Othwerwise; and array {v} must have
    {xmax-xmin+1} elements. */

void float_image_set_sample_row
  ( float_image_t *A, 
    int32_t c, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t y, 
    float v[] 
  );
  /* Stores the sample value {v[x-xmin]} into channel {c}, column {x}, row {y} of
    image {A}, for {x} in {xmin..xmax}.  
    
    If {xmin>xmax)} no samples are copied. Otherwise the range
    {xmin..xmax} must be contained in {0..NX-1}, where {NX} is the number
    of columns of {A}; and {v} must have {xmax-xmin+1} elements. */
    
void float_image_fill_sample_row 
  ( float_image_t *A,  
    int32_t c,  
    int32_t xmin,  
    int32_t xmax,  
    int32_t y,  
    float v 
  );
  /* Stores the sample value {v} into channel {c}, column {x}, row {y}
    of image {A}, for {x} in {xmin..xmax}.  
    
    If {xmin>xmax} no samples are set. Otherwise the range {xmin..xmax}
    must be contained in {0..NX-1}, where {NX} is the number of columns
    of {A}. */

/* ACESSING WHOLE ROWS OF PIXELS */

void float_image_get_pixel_row(float_image_t *A, int32_t xmin, int32_t xmax, int32_t y, float v[]);
  /* Sets {v[c + (x-xmin)*NC]} to the sample of channel {c}, column {x},
    row {y} of image {A}, for {c} in {0..NC-1} and {x} in {xmin..xmax};
    or to {NAN} if{ x} is not in {0..NX-1} -- where {NC} and {NX} are
    the number of channels and number of columns of {A}.
    
    If {xmin>xmax)} no samples are copied. Otherwise the range {xmin..xmax}
    must be contained in {0..NX-1}, where {NX} is the number of columns
    of {A}; and {v} must have {NC*(xmax-xmin+1)} elements. */

void float_image_set_pixel_row(float_image_t *A, int32_t xmin, int32_t xmax, int32_t y, float v[]);
  /* Stores the sample value {v[c + x*NC]} into channel {c}, column
    {x}, row {y} of image {A}, for {c} in {0..NC-1} and {x} in {xmin..xmax},
    where {NC} is the number of channels of {A}.
    
    If {xmin>xmax)} no samples are copied. Otherwise the range {xmin..xmax}
    must be contained in {0..NX-1}, where {NX} is the number of columns
    of {A}; and {v} must have {NC*(xmax-xmin+1)} elements. */
    
void float_image_fill_pixel_row(float_image_t *A, int32_t xmin, int32_t xmax, int32_t y, float v);
  /* Stores the sample value {v} into channel {c}, column {x}, row {y}
    of image {A}, for {c} in {0..NC-1} and {x} in {xmin..xmax},
    where {NC} is the number of channels of {A}.
    
    If {xmin>xmax)} no samples are set. Otherwise the range {xmin..xmax}
    must be contained in {0..NX-1}, where {NX} is the number of columns
    of {A}. */

/* ACCESSING WHOLE CHANNELS */

void float_image_fill_channel(float_image_t *A, int32_t c, float v);
  /* Stores {v} into all samples in channel {c} of image {A}. */

void float_image_set_channel(float_image_t *A, int32_t cA, float_image_t *V, int32_t cV);
  /* Copies all samples of channel {cV} of image {V}
    into channel {cA} of image {A}. */

void float_image_mix_channels
  ( double sA,
    float_image_t *A,
    int32_t cA, 
    double sB,
    float_image_t *B,
    int32_t cB, 
    float_image_t *R,
    int32_t cR
  );
  /* Sets each sample of channel {cR} of image {R} to {sA*vA + sB*vB},
    where {vA,vB} are the corresponding samples of channel {cA} of
    image {A} and channel {cB} of image {B}. It is safe to use {A==R}
    and/or {B==R}. */

/* ACCESSING THE WHOLE IMAGE */

void float_image_fill(float_image_t *A, float v);
  /* Stores {v} into all samples of image {A}. */
 
void float_image_fill_pixels(float_image_t *A, float v[]);
  /* Stores {v[c]} into all samples of channel {c} of image {A},
    for {c} in {0..NC-1}, where {NC} is the number of channels of {A}. */

void float_image_assign(float_image_t *A, float_image_t *V);
  /* Copies all samples of image {V} into image {A}.
    They must have the same dimensions and channel counts. */

/* ACCESSING RECTANGULAR SUB-IMAGES: */

void float_image_fill_rectangle
  ( float_image_t *A, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t ymin, 
    int32_t ymax, 
    float v
  );
  /* Stores {v} into all the samples of image {A} 
    in columns {xmin..xmax} and rows {ymin..ymax}. */ 

void float_image_fill_channel_rectangle
  ( float_image_t *A, 
    int32_t c, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t ymin, 
    int32_t ymax,
    float v
  );
  /* Stores {v} into all the samples of channel {c} of image {A} 
    in columns {xmin..xmax} and rows {ymin..ymax}. */ 
 
void float_image_fill_rectangle_pixels
  ( float_image_t *A, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t ymin, 
    int32_t ymax, 
    float v[]
  );
  /* Stores {v[0..NC-1]} into all pixels of image {A} 
    in columns {xmin..xmax} and rows {ymin..ymax}. */ 

void float_image_assign_channel_rectangle
  ( float_image_t *A, 
    int32_t cA, 
    int32_t xminA, int32_t xmaxA,  
    int32_t yminA, int32_t ymaxA,  
    float_image_t *V,  
    int32_t cV,  
    int32_t xminV, 
    int32_t yminV 
  );
  /* Copies into a rectangular region of channel {cA} of image {A} the
    corresponding samples from a rectangular region in channel {cV} of image {V}.
    The rectangle in {A} spans columns {xminA..xmaxA} and rows {yminA..ymaxA}.
    The rectangle in {V} has the same size but starts at column {xminV}
    and row {yminV}. All those samples must exist in the respective images. */ 

void float_image_assign_rectangle
  ( float_image_t *A,  
    int32_t xminA, int32_t xmaxA,  
    int32_t yminA, int32_t ymaxA,  
    float_image_t *V,  
    int32_t xminV,  
    int32_t yminV 
  );
  /* Copies into a rectangular region of image {A} the
    corresponding pixels from a rectangular region of image {V}.
    The rectangle in {A} spans columns {xminA..xmaxA} and rows {yminA..ymaxA}.
    The rectangle in {V} has the same size but starts at column {xminV}
    and row {yminV}. The images must have the same channel counts.
    All those pixels must exist in the respective images. */ 

/* PIXEL-WISE TRANSFORMATIONS */
    
void float_image_make_grayscale(float_image_t *A);
  /* Converts the image {A} to grayscale.  If {A} has a single channel,
    does nothing. If {A} has three channels, assumes that they are 
    the red, green and blue channels of the RGB model with linear
    enconding, and replaces all samples of each pixel by its brightness,
    computed with {frgb_get_Y}. If the number of channels is not 1 or 3,
    the procedure fails. */ 

void float_image_apply_gamma(float_image_t *A, int32_t c, double expo, double bias);
  /* Applies to channel {c} of image {A} the power correction with exponent {expo}
    and offset {bias}. See {sample_conv_gamma} for details. */ 

void float_image_log_scale(float_image_t *A, int32_t c, double bias, double vRef, double base);
  /* Converts every sample {v} in channel {c} of image {A} from linear
    to logarithmic scale, relative to value {vRef} and the given {base},
    by calling {sample_conv_log(v,bias,vRef,log(base))} (q.v.) for every
    sample {v} in channel {c} of {A}. */ 

void float_image_undo_log_scale(float_image_t *A, int32_t c, double bias, double vRef, double base);
  /* Undoes the effect of {float_image_log_scale}, by calling
    {sample_conv_undo_log(v, bias, vRef, log(base))} (q.v.) for every sample
    {v} in channel {c} of {A}. */ 

void float_image_rescale_samples(float_image_t *A, int32_t c, float a0, float a1, float z0, float z1);
  /* Applies to every sample of channel {c} of {A} the affine
    rescaling that takes the interval {[a0 _ a1]} to {[z0 _ z1]}. The
    samples are not clipped to either interval.
    
    The rescaling is defined as {v = z0 + (v-a0)*((z1-z0)/(a1-a0))}
    where {v} is the sample's value.
    
    There is no requirement that {a0<a1} or {z0<z1}. However, if {a0==a1}, 
    the rescaling may produce samples that are {±INF} or {NAN}. */

void float_image_square_samples(float_image_t *A, int32_t c);
  /* Replaces every sample in channel {c} of {A} by its square.
    In particular, leaves {+INF} and {NAN} alone, and maps {-INF} to {+INF}.
    Fails is {c} is not a valid channel index. */

double float_image_compute_sample_sum(float_image_t *A, int32_t c, int32_t *NS_P);
  /* Computes the sum of all samples in channel {c} of {A}.
    Ignores samples that are {±INF} or {NAN}. If {c} is not a valid
    channel index, returns 0.  If {NS_P} is not {NULL}, returns 
    in {*NS_P} the count of samples that were considered. */

double float_image_compute_squared_sample_sum(float_image_t *A, int32_t c, double avg, int32_t *NS_P);
  /* Computes the sum of {(s - avg)^2} over all samples {s} in channel {c} of {A}.
    Ignores samples that are {±INF} or {NAN}. If {c} is not a valid
    channel index, returns 0. If {NS_P} is not {NULL}, returns 
    in {*NS_P} the count of samples that were considered. */
    
void float_image_replace_nan_samples(float_image_t *A, int32_t c, float v);
  /* Replaces all {NAN} samples in channel {c} of {A} by {v}. */
    
/* PIXEL STATISTICS */

void float_image_compute_sample_avg_dev(float_image_t *A, int32_t c, double *avg, double *dev);
  /* Computes the average {*avg} and standard deviation {*dev} of all
    samples in channel {c} of image {A}. 
    
    Ignores samples that are {±INF} or {NAN}. The deviation uses the
    unbiased variance estimator formula (with {N-1} rather than {N} in
    the denominator). If {c} is not a valid channel index, assumes
    that all samples in that channel are all zero. If there are no
    valid samples, both parameters are set to zero. If there is only
    one valid sample, the average is that sample, and the deviation is
    set to zero. */

void float_image_update_sample_range(float_image_t *A, int32_t c, float *vMin, float *vMax);
  /* Expands the range {[*vMin _ *vMax]} so that it contains all
    finite sample values in channel {c} of image {A}.  Note that 
    {vMin,vMax} MUST be initialized by the client.
    
    Ignores samples that are {±INF} or {NAN}. If {c} is not a valid
    channel index, assumes that all samples in that channel are all
    zero. If there are no valid samples, does not modify {*vMax}
    and {*vMin}.  Either pointer may be null. */
    
float float_image_spectrum_max_sample(float_image_t *A, int32_t c, bool_t centered);
  /* Returns the maximum sample value in channel {c} the image {A}, which is assumed to 
    be a power spectrum image.  Ignores samples that are {±INF} or {NAN},
    and the sample corresponding to the zero-frequency (constant) Fourier 
    term. 
    
    The entries of {A} must be non-negative. 
    The zero-frequency term is assumed to be pixel {0,0} if {centered}
    is false, and {NX/2,NY/2} if {centered} is true.  Returns 0.0
    if there are no valid pixels in channel {c}, or {c} is an
    invalid channel index. */

/* IN-PLACE PIXEL REARRANGEMENT */

void float_image_flip_x(float_image_t *A, int32_t c, int32_t ix);
  /* Flips the contents of channel {c} of image {A} left-to-right, and
    cyclically shifts it in the X direction so that the original
    contents of column 0 ends up in column {ix}. Fails if {c} is an
    invalid channel index.
    
    The index {ix} and all column indices are implicitly reduced
    modulo the image width {nx}, to the range {0..nx-1}. Thus,
    if {nx == 7} and {ix == 9}, a row that originally contained the
    values {0 1 2 3 4 5 6} will contain {2 1 0 6 5 4 3}.  */

void float_image_flip_y(float_image_t *A, int32_t c, int32_t ix);
  /* Flips the contents of channel {c} of image {A} top-to-bottom, and
    cyclically shifts it in the Y direction so that the original
    contents of row 0 ends up in row {iy}. The effect is analogous to
    {float_image_flip_x} with swapped coordinate axes. Fails if {c} is
    an invalid channel index. */

void float_image_shift(float_image_t *A, int32_t c, int32_t dx, int32_t dy);
  /* Shifts the contents of channel {c} of image {A} by {dx} columns and 
    {dy} rows, cyclically.  
    
    Namely, the sample value that was originally stored in column {x}
    and row {y} will be moved to column {(x+dx) mod nx} and row
    {(y+dy) mod ny}, for all {x} and {y}; where {nx,ny} are the image
    dimensions, and {mod} means the mathematical remainder (always
    non-negative). */

/* GET/CHECK SIZE */
  
void float_image_get_size(float_image_t *A, int32_t *NC_P, int32_t *NX_P, int32_t *NY_P);
  /* Stores the number of channels, columns, and rows of {A} in
    {*NC_P}, {*NX_P}, and {*NY_P}, respectively. Addresses that are {NULL},
    are not set. */
  
void float_image_check_size(float_image_t *A, int32_t NC, int32_t NX, int32_t NY);
  /* Bombs out if {A} does not have exactly {NC} channels, {NX}
    columns, and {NY} rows.  If any of {NC,NX,NY} is negative,
    the corresponding parameter of {A} is not checked. */

/* INPUT/OUTPUT */

void float_image_write(FILE *wr, float_image_t *A);
  /* Writes the image {A} to {wr}, in ASCII floating point format.
    The file will start with a standard header line 
    "begin float_image_t (version {DATE})", 
    followed by three lines with "NC = {NUM}", "NX = {NUM}", "NY =
    {NUM}". Then follow the pixels, one per line with "{X} {Y}
    {SAMPLE}...". A blank line is inserted between two consecutive
    rows. The file ends with a footer "end float_image_t". */

float_image_t *float_image_read(FILE *rd);
  /* Reads an image from {rd}, in ASCII floating point format
    generated by {float_image_write}. */

/* DEBUGGING */

void float_image_debug_pixel(char *label, double x, double y, int32_t chns, float f[], char *tail);
  /* Prints the point {p} and the float-valued samples
    {f[0..chns-1]} to {stderr}, tagged with the string {label} and
    terminated by the string {tail}. */

void float_image_debug_double_pixel(char *label, double x, double y, int32_t chns, double v[], char *tail);
  /* Prints the point {(x,y)} and the float-valued samples
    {v[0..chns-1]} to {stderr}, tagged with the string {label} and
    terminated by the string {tail}. */

#endif
