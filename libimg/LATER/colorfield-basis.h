/* Representation, evaluation and fitting of color fields. */
/* Last edited on 2024-12-04 23:19:20 by stolfi */

#ifndef colorfield_H
#define colorfield_H

#include <stdio.h>
#include <stdint.h>

#include <FBas.h>

#define MAXCHANNELS 3
  
typedef struct Triplet { double c[MAXCHANNELS]; } Triplet; 

#define Black   (Triplet){{0.0, 0.0, 0.0}}
#define White   (Triplet){{1.0, 1.0, 1.0}}
#define Ones    (Triplet){{1.0, 1.0, 1.0}}
#define NoColor (Triplet){{NAN, NAN, NAN}}

int32_t eq_triplet(Triplet a, Triplet b, int32_t chns);
  /* TRUE iff {a.c[i]} equals {b.c[i]} for {i} in {0..chns-1}. */

void print_triplet(FILE *f, char *pref, int32_t chns, Triplet *fv, char *suff);
  /* Print the given triplet to file {f}, surrounded by the given {pref} and
    {suff} strings. */

typedef struct ColorField
  { FBas *bas;      /* The approximating basis. */
    int32_t cols;        /* Number of columns in image. */
    int32_t rows;        /* Number of rows in image. */
    int32_t nt;          /* Number of relevant elements in the field. */
    int32_t *index;      /* Indices of significant basis elements. */
    Triplet *wt;     /* Color coefficients of those elements. */
    int32_t additive;    /* TRUE to add, FALSE to multiply. */
  } ColorField;
  /* Describes a color field {f}, i.e. a function from pixel indices
    {col,row} to colors. The field is obtained by linear combination of
    a fixed fixed list of functions {bas[0..nb-1]} (the /basis/).
    The basis is determined by its {bas} and its size {nb}.
  
    If {additive = TRUE}, each component of {f(p)} at a pixel {p} is
    the sum {wfun(p)} of {nt} basis elements {bas[index[0..nt-1]](p)},
    with weights {w[0..nt-1]}.
    
    If {additive = FALSE}, the field value is {f(p) = exp(wfun(p))},
    where {wfun(p)} is defined as above.
    
    The basis elements are evaluated with relative coordinates whose
    range is {[-1 _ +1]}, where {(-1,-1)} is the pixel at the top left
    corner of the image, and {(+1,+1)} at the bottom right*/ 
  
int32_t field_is_constant(ColorField *cf, Triplet color, int32_t chns);
  /* TRUE if the field {cf} is constant, and its 
    first {chns} channels match {color[0..chns-1]}.
    WARNING: Assumes that the basis element with index 0
    is the constant function 1. */

Triplet eval_color_field(ColorField *cf, int32_t col, int32_t row, int32_t chns);
  /* Evaluates the color field {cf} at the pixel in column {col} and
    row {row} of the image. Only channels {[0..chns-1]} are
    evaluated. */

/* COLORFIELD DATA SAMPLES */
    
typedef struct PixPos { int32_t c[2]; } PixPos; 

typedef struct ColorData  /* Color field specifications from command line. */
  { int32_t np;           /* Number of sample points. */
    PixPos *pt;       /* Pixel indices of sample points. */
    Triplet *color;   /* Corresponding colors. */
  } ColorData;
  /* A {ColorData} record contains user-oriented specification of a
    color field. The field is defined by a set of sample points
    {pt[i]} and the corresponding color values {color[i]}. */

/* FITTING COLOR FIELDS TO DATA SAMPLES */
    
ColorField *compute_black_field
  ( ColorData *cd,    /* User specs of color field. */
    FBas *bas,        /* The approximation basis. */ 
    int32_t chns,         /* Color channels to use. */
    int32_t cols,         /* Number of columns in image. */
    int32_t rows,         /* Number of rows in image. */
    double tol,       /* Discard terms that are less than this. */
    Triplet *inGamma  /* Gamma of color samples, per channel. */
  );
  /* Computes an additive gamma-corrected black reference
    {*blackField}, from the space generated by {bas}, that fits the
    data {cd}. The sampled colors in {cd} will be interpreted
    according to {inGamma}. */

ColorField *compute_white_field
  ( ColorData *cd,     /* User specs of color field. */
    FBas *bas,         /* The approximation basis. */ 
    int32_t chns,          /* Color channels to use. */
    int32_t cols,          /* Number of columns in image. */
    int32_t rows,          /* Rows in image. */
    double tol,        /* Discard terms that are less than this. */
    Triplet *inGamma,  /* Gamma of color samples, per channel. */
    ColorField *black, /* A previously computed {black} color field. */
    double minw        /* Minimum intensity for {white-black}. */
  );
  /* Computes a gamma- and black-corrected multiplicative white
    reference field, from the space generated by {bas}, that fits the
    data {cd}. The color values in {cd} will be interpreted according
    to {inGamma}), and the previously computed {black} field will be
    subtracted from them. The result is increased if needed to be at
    least {minw} in each channel. */

ColorField *fit_color_field
  ( FBas *bas,      /* Basis to use. */  
    int32_t np,         /* Number of samples. */
    PixPos *pt,     /* Sample positions. */ 
    Triplet *color, /* Sample values. */  
    int32_t chns,       /* Color channels to use. */
    int32_t cols,       /* Columns in image. */
    int32_t rows,       /* Rows in image. */
    double tol      /* Discard terms that are less than this. */
  );
  /* Computes an additive color field of the given {bas} that best matches 
    a set of samples {pt}. Assumes that the colors {color} have already been
    corrected for gamma and other applicable factors.  Discards terms 
    whose absolute value is less than {tol}. */

/* DEBUGGING AIDS */

extern int32_t cfld_debug;    /* Set to TRUE to print debugging info */

void cfld_debug_triplet(char *label, int32_t col, int32_t row, int32_t chns, Triplet *fv, char *tail);
  /* Print the given triplet if the global flag {debug} is true. */

#endif
