#ifndef image_window_op_H
#define image_window_op_H

/* {image_window_op.h} - local neighborhood image operators. */
/* Created on 2012-10-25 by J. Stolfi, UNICAMP */
/* Last edited on 2024-12-05 10:30:38 by stolfi */

#include <stdint.h>

#include <bool.h>

typedef enum
  { image_window_op_AVERAGE,    /* Weighted average. */
    image_window_op_VARIANCE,   /* Weighted variance from average. */
    image_window_op_DEVIATION,  /* Weighted deviation from average. */
    /* Differential ops based on quadratic poly fit: */
    image_window_op_F,          /* Function value. */
    image_window_op_FX,         /* Horizontal derivative. */
    image_window_op_FY,         /* Vertical derivative. */
    image_window_op_FXX,        /* Second derivative wrt {x}. */
    image_window_op_FXY,        /* Mixed second derivative {x,y}. */
    image_window_op_FYY,        /* Second derivative wrt {y}. */
    image_window_op_FXXY,       /* (Pseudo) mixed derivative wrt {x,x,y}. */
    image_window_op_FXYY,       /* (Pseudo) mixed derivative wrt {x,y,y}. */
    image_window_op_FXXYY,      /* (Pseudo) mixed derivative wrt {x,x,y,y}. */
    image_window_op_GRADIENT,   /* Gradient modulus. */
    image_window_op_LAPLACIAN,  /* Laplacian {fxx+fyy}. */
    image_window_op_ORTHICITY,  /* Orthicity. */
    image_window_op_ELONGATION, /* Elongation. */
    /* Operators based on weighted least squares degree 1 poly fit: */    
    image_window_op_LINF,       /* Constant term of LSQ affine fit. */
    image_window_op_LINFX,      /* Coefficient of {x} in LSQ affine fit. */
    image_window_op_LINFY,      /* Coefficient of {x} in LSQ affine fit. */
    image_window_op_LINVAR,     /* Weighted variance of error from LSQ affine fit. */
    image_window_op_LINDEV,     /* Weighted deviation of error from LSQ affine fit. */
    /* Sentinel: */    
    image_window_op_NUM         /* Number of valid kinds. Not itself a valid kind. */
  } image_window_op_t;
  /* Kind of operator. */
  
#define image_window_op_HELP_OPERATORS \
  "average | variance | deviation \\\n" \
  "    f | fx | fy | fxx | fxy | fyy | fxxy | fxyy | fxxyy | \\\n" \
  "    gradient | laplacian | orthicity | elongation \\\n" \
  "    linf | linfx | linfy | linvar | lindev"

#define image_window_op_INFO_OPERATORS \
  "Average and deviation\n" \
  "\n" \
  "  The {average} operator is the weighted average of the window pixels, with weight 4 for the" \
  " center pixel, 2 for each of its four nearest neighbors, and 1 for the corner" \
  " pixels.  Namely\n" \
  "\n" \
  "      {average}    {((4*Woo + 2*(Wpo+Wmo+Wop+Wom) + (Wpp+Wpm+Wmp+Wmm))/16}\n" \
  "\n" \
  "  The {variance} operator is the weighted average value of {(W[x+r,y+s]-avg)^2}, the" \
  " square of the differences between each sample value and that {average}, using the" \
  " same 1-2-4 weights.  The {deviation} operator is the square root of" \
  " {variance}, namely the weighted root mean square displacement of samples from the {average}.\n" \
  "\n" \
  "Differential operators\n" \
  "\n" \
  "  The \"differential\" operators {f,fx,fy,fxx,fxy,fyy,fxxy,fxyy,fxxyy} are" \
  " conceptually defined by  fitting a polynomial of maximum" \
  " (tensor) degree 2 to the sample values in the window, and" \
  " differentiating that polynomial.  Namely, the procedure" \
  " determines nine coefficientst {A,B,...H,I} such that the" \
  " sample value {W[x,y]} is the value of the polynomial\n" \
  "\n" \
  "  {P(x,y) = A + B*x + C*y + D*x^2 + E*x*y + F*x^2 + G*x^2*y + H*x*y^2 + I*x^2*y^2}\n" \
  "\n" \
  "  where {x} is the column index and {y} is the row index," \
  " both in {-1.0,+1}. Then {f,fx,fy,fxy} are {A,B,C,E}, {fxx,fyy,fxxy,fxyy} are" \
  " twice {D,F,G,H}, and {fxxyy} is four times {I}.  This model" \
  " implies that other third and fourth" \
  " derivatives {fxxx,fyyy,fxxxx,fxxxy,fxyyy,fyyyy} are" \
  " all zero.  In terms of the window samples, the formulas for these" \
  " derivatives turn out to be:\n" \
  "\n" \
  "      {f}          {Woo}\n" \
  "      {fx}         {(Wpo-Wmo)/2}\n" \
  "      {fy}         {(Wop-Wom)/2}\n" \
  "      {fxx}        {Wmo+Wpo-2*Woo}\n" \
  "      {fxy}        {(Wpp-Wmp-Wpm+Wmm)/4}\n" \
  "      {fyy}        {Wom+Wop-2*Woo}\n" \
  "      {fxxy}       {(Wpp-2*Wop+Wmp-Wpm+2*Wom-Wmm)/2}\n" \
  "      {fxyy}       {(Wpp-2*Wpo+Wpm-Wmp+2*Wmo+Wmm)/2}\n" \
  "      {fxxyy}      {Wpp+Wpm+Wmp+Wmm-2*Wpo-2*Wmo-2*Wop-2*Wom+4*Woo}\n" \
  "\n" \
  "  where {Woo = W[x,y]}, {Wpo = W[x+1,y]}, {Wmo = W[x-1,y]}," \
  " {Wop = W[x,y+1]}, {Wom = W[x,y-1]}, {Wpp = W[x+1,y+1]}," \
  " {Wmp = W[x-1,y+1]}, etc.\n" \
  "\n" \
  "  Note that the {y} axis is assmumed to point towards increasing row" \
  " indices, which is 'down' if the image is displayed in the usual" \
  " way (with pixel [0,0] at top left).  If the {y} is assumed to point 'up', then" \
  " the values of {fy}, {fxy}, and {fxxy} should be negated.\n" \
  "\n" \
  "  Note also that the formulas" \
  " will give incorrect results if the image samples are affected" \
  " by random noise, or come from a polynomial or series with" \
  " terms with degree 3 or more in {x} or in {y}.\n" \
  "\n" \
  "Additional diferential operators\n" \
  "  Other \"differential\" operators are obtained by combining those nine in various ways:\n" \
  "\n" \
  "      {gradient}   {hypot(fx,fy)}\n" \
  "      {laplacian}  {fxx+fyy}\n" \
  "      {orthicity}  {fxx-fyy}\n" \
  "      {elongation} {hypot(2*fxy,fxx-fyy)}\n" \
  "\n" \
  "\n" \
  "  The formulas for {laplacian} and {orthicity} turn out to be:\n" \
  "\n" \
  "      {laplacian}  {Wop+Wom+Wmo+Wpo - 4*Woo}\n" \
  "      {orthicity}  {Wmo+Wpo-Wom-Wop}\n" \
  "\n" \
  "  The names \"orthicity\" and \"elongation\" are made up.  In the" \
  " model of images as continuous functions from {R^2} to {R}, with" \
  " the derivatives of Calculus, the" \
  " {laplacian} and {elongation} are the two second-order" \
  " local operators that are invariant under image rotation around the point.\n" \
  "\n" \
  "Smoothed differential operators\n" \
  "\n" \
  "  There is also a \"smoothed\" version of the differential" \
  " operators, that instead fits the following polynomial to" \
  " the window values:\n" \
  "\n" \
  "    {Q(x,y) = A + B*x + C*y + D*(2*x^2-1) + E*x*y + F*(2*y^2-1) + \n" \
  "              G*(2*x^2-y^2)*y + H*(2*y^2-x^2)*x + I*(1 - 2*(x^2 - y^2)^2) \n" \
  "\n" \
  "  and then computes the derivatives of this polynomial.  Namely, {f} is" \
  " {A-D-F+I}, {fx,fy,fxy} are {B,C,E}, {fxx,fyy,fxxy,fxyy} are" \
  " four times {D,F,G,H}, and {fxxyy} is 16 times {I}.  In terms" \
  " of the window samples, the formulas turn out to be:\n" \
  "\n" \
  "      smoothed {fx}    {((Wpp+2*Wpo+Wpm) - (Wmp+2*Wmo+Wmm))/8}\n" \
  "      smoothed {fy}    {((Wpp+2*Wop+Wmp) - (Wpm+2*Wom+Wmm))/8}\n" \
  "      smoothed {fxx}   {((Wmp+2*Wmo+Wmm) + (Wpp+2*Wpo+Wpm) - 2*(Wop+2*Woo+Wom))/4}\n" \
  "      smoothed {fyy}   {((Wpm+2*Wom+Wmm) + (Wpp+2*Wop+Wmp) - 2*(Wpo+2*Woo+Wmo))/4}\n" \
  "\n" \
  "  The smoothed version of {f} turns out to be the same as the non-smoothed" \
  " one, namely the center sample.  The smoothed versions of {fx} and {fy} are" \
  " the components of the Sobel gradient operator.  The smoothed versions" \
  " of {fxy,ffxxy,fxyy,fxxyy} turn out to be the same as the non-smooth ones.\n" \
  "\n" \
  "  Unlike the {P} polynomial, {Q} has additional non-zero third derivatives at the" \
  " origin, namely {fxxx = -6*H = -(3/2)*fxyy} and" \
  " {fyyy = -6*G = -(3/2)*fxxy}.\n" \
  "\n" \
  "   The smoothed versions of {gradient}, {laplacian}, {orthicity}, and {elongation} are" \
  " computed by the same formulas above, but using the" \
  " smoothed versions of {fx},{fy}, {fxx}, and {fyy}.  The smoothed " \
  " versions of {orthicity}, and {elongation} turn out to" \
  " be identical to the unsmoothed versions.   The Laplacian" \
  " formula {fxx+fyy}, computed with smoothed derivatives, comes out to\n" \
  "\n" \
  "      {fxx+fyy = (Wpp+Wmm+Wmp+Wmm - 4*Woo)/2}\n" \
  "\n" \
  "Linear fit operators\n" \
  "\n" \
  "  Another set of window operators approximates the window samples by a" \
  " first-degree polynomial {A + B*x + C*y}, using weighted least-squares fitting" \
  " with the 1-2-4 weights of {average}.  The derivatives of this polynomial -- that" \
  " is, the coefficients {A,B,C} -- are the the operators {linf,linfx,linfy}.  The" \
  " operator {lingrad} is {hypot(linfx,linfy).  The operator {linvar} is defined" \
  " as the weighted average of {W[x,y]-(A + B*x + C*Y)} {lindev} is the root" \
  " mean square devation of the window samples from that polynomial\n" \
  "\n" \
  "  The {linf} coefficient turns out to be" \
  " the same as {average}, and {linfx,linfy} turn out to be the same" \
  " as the smoothed first-derivative operators {fx,fy}.\n" \
  "\n" \

void image_window_op_get_range(image_window_op_t op, bool_t smoothed, bool_t squared, double *loP, double *hiP);
  /* Sets {*loP} and {*hiP} to the range of the operator {op} with the given {smoothed} and {squared} options,
    assuming the samples in the window are all in {[0_1]}. */

#define image_window_op_INFO_OP_RANGES \
  "    The natural ranges (before rescaling) for the operators, normal" \
  " and squared, are:\n" \
  "\n" \
  "      Operator      Plain                Squared\n" \
  "\n" \
  "      {average}     {[0 _ 1]}            {[0 _ 1]}\n" \
  "      {variance}    {[0 _ 1/4]}          {[0 _ 1/16]}\n" \
  "      {deviation}   {[0 _ 1/2]}          {[0 _ 1/4]}\n" \
  "      {f}           {[0 _ 1]}            {[0 _ 1]}\n" \
  "      {fx,fy}       {[-1/2 _ +1/2]}      {[0 _ 1/4]}\n" \
  "      {fxx,fyy}     {[-2 _ +2]}          {[0 _ 4]}\n" \
  "      {fxy}         {[-1/2 _ +1/2]}      {[0 _ 1/2]}\n" \
  "      {fxxy,fxyy}   {[-2 _ +2]}          {[0 _ 4]}\n" \
  "      {fxxyy}       {[-8 _ 8]}           {[0 _ 64]}\n" \
  "      {gradient}    {[0 _ +sqrt(1/2)]}   {[0 _ 1/2]}\n" \
  "      {laplacian}   {[-4 _ +4]}          {[0 _ 16]}\n" \
  "      {orthicity}   {[-2 _ +2]}          {[0 _ 4]}\n" \
  "      {elongation}  {[0 _ sqrt(2)]}      {[0 _ 2]}\n" \
  "      {linf}        {[0 _ 1]}            {[0 _ 1]}\n" \
  "      {linfx,linfy} {[-1/2 _ +1/2]}      {[0 _ 1/4]}\n" \
  "      {linvar}      {[0 _ 1/4]}          {[0 _ 1/16]}\n" \
  "      {lindev}      {[0 _ 1/2]}          {[0 _ 1/4]}\n" \
  "\n" \
  "    The smoothed versions of the operators have the same" \
  " ranges, except\n" \
  "\n" \
  "      Operator             Plain                Squared\n" \
  "\n" \
  "      smoothed {gradient}  {[0 _ sqrt(5)/4]}    {[0 _ 0.3125]}\n" \
  "      smoothed {laplacian} {[-2 _ +2]}          {[0 _ 4]}\n"

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
    
/* INDIVIDUAL OPERATORS

  These procedures apply the relevant operator to the array {smp} of
  window samples and return the output sample value for that pixel. The
  meaning of the parameters {nx} and {ictr} are as in
  {image_window_op_apply}. */

double image_window_op_average(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_variance(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_deviation(int32_t ictr, int32_t nx, double smp[]);
  /* These procedures compute the weighted {average}, {variance}, and {deviation} of the
    samples in the window.  There are no smoothed versions. */

double image_window_op_f(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fx_unsmoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fy_unsmoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fxx_unsmoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fxy(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fyy_unsmoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fxxy(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fxyy(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fxxyy(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_gradient_unsmoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_gradient_unsmoothed_squared(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_laplacian_unsmoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_orthicity(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_elongation(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_elongation_squared(int32_t ictr, int32_t nx, double smp[]);
  /* The differential operators derived from a quadratic polynomial fit. */

double image_window_op_fx_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fy_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fxx_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_fyy_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_gradient_smoothed(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_gradient_smoothed_squared(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_laplacian_smoothed(int32_t ictr, int32_t nx, double smp[]);
  /* The smoothed differential operators that are different from the unsmoothed versions. */

double image_window_op_linvar(int32_t ictr, int32_t nx, double smp[]);
double image_window_op_lindev(int32_t ictr, int32_t nx, double smp[]);
  /* The operators based on weighted least squares fit of an affine function.
    For {linf}, use {average}. For {linfx} and {lnfy}, use
    the smoothed {fx} and {fy}. */

/* OPERATOR NAMES */

const char *image_window_op_to_string(image_window_op_t op);
  /* Returns a string with the name of the operator {op} in lower case, 
     namely "f", "fx", "laplacian", "average", etc.. */

image_window_op_t image_window_op_from_string(const char *chop);
  /* Converts an operator name string {chop} in lower case, 
     namely "f", "fx", "laplacian", "average", etc.,
     to the corresponding operator. */

#endif
