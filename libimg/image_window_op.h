#ifndef image_window_op_H
#define image_window_op_H

/* {image_window_op.h} - local neighborhood image operators. */
/* Created on 2012-10-25 by J. Stolfi, UNICAMP */
/* Last edited on 2020-11-29 10:01:25 by jstolfi */

#define _GNU_SOURCE_
#include <stdint.h>

#include <bool.h>

typedef enum
  { image_window_op_IDENT,      /* Identity operator. */
    image_window_op_DX,         /* Horizontal derivative. */
    image_window_op_DY,         /* Vertical derivative. */
    image_window_op_DXX,        /* Second derivative wrt {x}. */
    image_window_op_DXY,        /* Mixed second derivative. */
    image_window_op_DYY,        /* Second derivative wrt {y}. */
    image_window_op_GRADIENT,   /* Gradient modulus. */
    image_window_op_LAPLACIAN,  /* Laplacian {dxx+dyy}. */
    image_window_op_ORTHICITY,  /* Orthicity {dxx-dyy}. */
    image_window_op_ELONGATION, /* Elongation modulus {hypot(2*dxy,dxx-dyy)}. */
    image_window_op_AVERAGE,    /* Weighted average. */
    image_window_op_DEVIATION,  /* Weighted deviation from average. */
    image_window_op_NUM_VALUES  /* Number of valid kinds. Not itself a valid kind. */
  } image_window_op_t;
  /* Kind of operator. */
  
#define image_window_op_HELP_OPERATORS \
  "ident | dx | dy | dxx | dxy | dyy | gradient | laplacian | orthicity | elongation | average | deviation"

#define image_window_op_INFO_OPERATORS \
  "    The \"ident\" operation is the identity operator, simply the central window" \
  " sample. The operators \"dx\", \"dy\", \"dxx\", \"dxy\", \"dyy\" are the spatial" \
  " derivatives of the input image.  Each derivative are computed" \
  " by the simplest unbiased finite difference\n" \
  " formula.  Specifically, the sample {V[x,y]}, on column {x} and row {y} of {V}, is" \
  " computed from the input image {U} as follows:\n" \
  "\n" \
  "      {ident}      {Uoo}\n" \
  "      {dx}         {(Upo-Umo)/2}\n" \
  "      {dy}         {(Uop-Uom)/2}\n" \
  "      {dxx}        {Umo+Upo-2*Uoo}\n" \
  "      {dyy}        {Uom+Uop-2*Uoo}\n" \
  "      {dxy}        {(Upp-Ump-Upm+Umm)/4}\n" \
  "      {gradient}   {hypot(dx,dy)}\n" \
  "      {laplacian}  {dxx+dyy = Uop+Uom+Umo+Upo - 4*Uoo}\n" \
  "      {orthicity}  {dxx-dyy = Umo+Upo-Uom-Uop}\n" \
  "      {elongation} {hypot(2*dxy,dxx-dxy)}\n" \
  "      {average}    {weighted average of 3x3 window}\n" \
  "      {deviation}  {root weighted mean square deviation from average}\n" \
  "\n" \
  "  where {Uoo = U[x,y]}, {Upo = U[x+1,y]}, {Umo = U[x-1,y]}, {Uop = U[x,y+1]}, {Uom = U[x,y-1]}, {Upp = U[x+1,y+1]}, {Ump = U[x-1,y+1]}, etc.\n" \
  "\n" \
  "    The names \"orthicity\" and \"elongation\" are made up.  The" \
  " Laplacian and elongation are the two second-order" \
  " local operators that (in the continuous model) are invariant under image rotation.  The" \
  " \"average\" uses weights 4 for the" \
  " center pixel, 2 for its four nearest neighbors, and 1 for the corner" \
  " pixels.  The variance is the weighted mean of {(U[x+r,y+s]-avg)^2}, using" \
  " the same weights; the deviation is the square root of the variance.\n" \
  "\n" \
  "    There is also a \"smoothed\" version of the differential operators that is applied to an implicitly version of the image smoothed by the 4-2-1 mask above.  Namely:\n" \
  "\n" \
  "      smoothed {ident} {((4*Uoo + 2*(Upo+Umo+Uop+Uom) + (Upp+Upm+Ump+Umm))/16}\n" \
  "      smoothed {dx}    {((Upp+2*Upo+Upm) - (Ump+2*Umo+Umm))/8}\n" \
  "      smoothed {dy}    {((Upp+2*Uop+Ump) - (Upm+2*Uom+Umm))/8}\n" \
  "      smoothed {dxx}   {((Ump+2*Umo+Umm) + (Upp+2*Upo+Upm) - 2*(Upo+2*Uoo+Umo))/4}\n" \
  "      smoothed {dyy}   {((Upm+2*Uom+Umm) + (Upp+2*Uop+Ump) - 2*(Uop+2*Uoo+Uom))/4}\n" \
  "      smoothed {dxy}   {(Upp-Ump-Upm+Umm)/4}\n" \
  "\n" \
  "    Note that the smoothed versions of \"dx\" and \"dy\" are the Sobel gradient components, and that the smoothed \"ident\" is the same as \"average\".  The smoothed versions of \"gradient\", \"laplacian\", \"orthicity\", and \"elongation\" are computed from the derivatives by the same formulas above, but using the smoothed versions of \"dx\",\"dy\", \"dxx\", and \"dxy\".   The smoothed versions of \"dxy\", \"orthicity\", and \"elongation\" turn out to be identical to the unsmoothed versions.   The smoothed Laplacian simplifies to\n" \
  "\n" \
  "      smoothed {laplacian}  {dxx+dyy = (Upp+Umm+Ump+Umm - 4*Uoo)/2}\n" \
  "\n" \
  "    The smoothed versions of \"average\" and \"deviation\" are identical to the unsmoothed ones, by definition."

#define image_window_op_INFO_OP_RANGES \
  "    The natural ranges (before rescaling) for the operators, normal and squared, are:\n" \
  "\n" \
  "      Operator     Plain                Squared\n" \
  "\n" \
  "      {ident}      {[0 _ 1]}            {[0 _ 1]}\n" \
  "      {dx}         {[-1/2 _ +1/2]}      {[0 _ 1/4]}\n" \
  "      {dy}         {[-1/2 _ +1/2]}      {[0 _ 1/4]}\n" \
  "      {dxx}        {[-2 _ +2]}          {[0 _ 4]}  \n" \
  "      {dyy}        {[-2 _ +2]}          {[0 _ 4]}  \n" \
  "      {dxy}        {[-1/2 _ +1/2]}      {[0 _ 1/2]}\n" \
  "      {gradient}   {[0 _ +sqrt(1/2)]}   {[0 _ 1/2]}\n" \
  "      {laplacian}  {[-4 _ +4]}          {[0 _ 16]} \n" \
  "      {orthicity}  {[-2 _ +2]}          {[0 _ 4]} \n" \
  "      {elongation} {[0 _ sqrt(5)]}      {[0 _ 5]}  \n" \
  "      {average}    {[0 _ 1]}            {[0 _ 1]}  \n" \
  "      {deviation}  {[0 _ 1/2]}          {[0 _ 1/4]}\n" \
  "\n" \
  "    The smoothed versions of the operators have the same ranges.\n" \

double image_window_op_apply
  ( image_window_op_t op, 
    bool_t smoothed, 
    bool_t squared, 
    int32_t ictr, 
    int32_t nx, 
    double smp[]
  );
  /* Applies the operator {op} to the window samples in the array {smp}, 
    assumed to be linearized by rows.  If {smoothed} is true, computes
    the smoothed version.  If {squared} is true, returns the square of the result.
    (For some operators, this option is faster than computing the unsquared
    operator and then squaring it).   
    
    The parameter {ictr} must be the index of the central window sample
    in the array {smp}, and {nx} must be the number of columns in
    that array. The array {smp} must have enough rows and columns to
    include all samples needed by the operator.
    
    If any sample {smp[k]} used by the operator is infinite, the result
    may be {±INF} or {NAN}. If any {smp[k]} used by the operator is
    {NAN}, the result is likely to be {NAN}. */

void image_window_op_get_window_size
  ( image_window_op_t op, 
    bool_t smoothed,
    int32_t *nxP, 
    int32_t *nyP, 
    int32_t *ictrP
  );
  /* Returns in {*nxP} and {*nyP} the minimum number of cols and
    rows needed to compute the operator {op}, with the given {smoothed}
    option. These numbers will always be odd. Also sets {*ictrP} to the
    index of the central sample in the linearized window. */

void image_window_op_get_range(image_window_op_t op, bool_t squared, double *loP, double *hiP);
  /* Sets {*loP} and {*hiP} to the range of the operator {op} with the given {squared} option,
    assuming inputs in {[0_1]}. */

double image_window_op_ident(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_dx(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_dy(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_dxx(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_dxy(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_dyy(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_gradient(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_gradient_squared(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_laplacian(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_orthicity(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_elongation(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_elongation_squared(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_average(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_deviation(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_deviation_squared(int32_t ictr, int32_t nx, double smp[]);
  /* These procedures apply the relevant operator to the
    array {smp} of window samples and return the output sample value for
    that pixel.  The parameters have the same meaning as in */

double image_window_op_dx_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_dy_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_dxx_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_dyy_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_gradient_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_gradient_smoothed_squared(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_laplacian_smoothed(int32_t ictr, int32_t nx, double smp[]);
  /* Smoothed operators that are different from the unsmoothed versions. */

const char *image_window_op_to_string(image_window_op_t op);
  /* Returns a string with the name of the operator {op} in lower case, 
     namely "ident", "dx", "laplacian", "average", etc.. */

image_window_op_t image_window_op_from_string(const char *chop);
  /* Converts an operator name string {chop} in lower case, 
     namely "ident", "dx", "laplacian", "average", etc.,
     to the corresponding operator. */

#endif
