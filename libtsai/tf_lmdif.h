/* Last edited on 2022-10-20 07:33:54 by stolfi */

#ifndef tf_lmdif_H
#define tf_lmdif_H

/* Solves or minimizes the sum of squares of m nonlinear
  functions of n variables.  From public domain Fortran version
  of Argonne National Laboratories MINPACK
  
  C translation by Steve Moshier

  Adapted by J.Stolfi, UNICAMP, 2015-10-04.
*/

#include <stdint.h>

typedef void tf_lmdif_fcn_t (int32_t m, int32_t n, double x[], double fvec[], int32_t *iflag);
  /* Type of the user-supplied procedure which
    calculates the functions for {lmdif}.  It receives an argument 
    vector {x[0..n-1]} and should store the nonlinear 
    functions in {fvec[0..m-1]}.  If {*iflag} is 0 on input,
    then the user has requsted printing (see {nprint} argument to 
    {lmdif} below).  The procedure may set {*iflag}
    to a negative integer to terminate {lmdif}. */

int32_t lmdif
  ( tf_lmdif_fcn_t *fcn,
    int32_t m,
    int32_t n,
    double x[],
    double fvec[],
    double ftol,
    double xtol,
    double gtol,
    int32_t maxfev,
    double epsfcn,
    double diag[],
    int32_t mode,
    double factor,
    int32_t nprint,
    int32_t *info,
    int32_t *nfev,
    double fjac[],
    int32_t ldfjac,
    int32_t ipvt[],
    double qtf[],
    double wa1[],
    double wa2[],
    double wa3[],
    double wa4[]
  );
  /* Minimizes the sum of the squares of {m} nonlinear functions in
    {n} variables by a modification of the Levenberg-Marquardt algorithm. 
    The user must provide a function {fcn} which calculates the functions.
    The Jacobian is then calculated by a forward-difference approximation.
   
      {fcn} is the procedure that evaluates the nonlinear functions
         whose sum of squares is to be minimized. The {fcn} procedure may
         set {*iflag} to a negative integer to terminate {lmdif}.
   
      {m} is the number of functions.
   
      {n} is the number of independent variables. It must not exceed {m}.
   
      {x} is an array of length {n}. on input {x[0..n-1]} must contain
          an initial estimate of the solution vector. on output {x[0..n-1]}
          contains the final estimate of the solution vector.
   
      {fvec} is an array of length {m}.  On output, {fvec[0..m-1]} will contain
          the functions evaluated at the output {x}.
   
      {ftol} is a non-negative value, the relative error desired
           in the sum of squares.  Termination
           occurs when both the actual and predicted relative
           reductions in the sum of squares are at most {ftol}.
   
      {xtol} is a nonnegative value, the relative error desired in the
           approximate solution.  Termination occurs when the relative
           error between two consecutive iterates is at most {xtol}.

      {gtol is a nonnegative value, the orthogonality
           desired between the function vector and the columns
           of the Jacobian.  Termination occurs when the cosine of
           the angle between {fvec} and any column of the jacobian
           is at most {gtol} in absolute value.
   
      {maxfev} is a positive integer.  Termination
          occurs when the number of calls to {fcn} is at least
          {maxfev} by the end of an internal iteration of {lmdif}.
   
      {epsfcn} is an input variable used in determining a suitable
          step length for the forward-difference approximation.  This
          approximation assumes that the relative errors in the
          functions are of the order of {epsfcn}. if {epsfcn} is less
          than the machine precision, it is assumed that the relative
          errors in the functions are of the order of the machine
          precision.
   
      {diag} is an array of length {n}.  If {mode = 1} (see
          below), {diag[0..n-1]} are set internally. 
          If {mode = 2}, {diag[0..n-1]} must contain positive
          entries that serve as multiplicative scale 
          factors for the variables.
   
      {mode} is an integer flag that specifies how the argument
          variables are scaled.  If {mode = 1}, the
          variables will be scaled internally.  if {mode = 2},
          the scaling is specified by the {diag}.  Other
          values of {mode} are equivalent to {mode = 1}.
   
      {factor} is a positive value, used in determining the
          initial step bound.  This bound is set to the product of
          {factor} and the Euclidean norm of {diag*x} if nonzero, or else
          to {factor} itself.  In most cases {factor} should lie in the
          interval {(0.1 _ 100)}; 100 is a generally recommended value.
   
      {nprint} is an integer input variable that enables controlled
          printing of iterates if it is positive.  In this case,
          {fcn} is called with {iflag = 0} at the beginning of the first
          iteration and every {nprint} iterations thereafter and
          immediately prior to return, with {x} and {fvec} available
          for printing. if nprint is not positive, no special calls
          of {fcn} with {iflag = 0} are made.
   
      {info} is the address of an integer output variable.  If the execution
          was terminated because a call to {fcn} set {*iflag} to a negative value,
          {info} is set to that (negative) value.  Otherwise,
          {info} is set as follows:
   
            {info = 0}  improper input parameters.
     
            {info = 1}  both actual and predicted relative reductions
                      in the sum of squares are at most {ftol}.
     
            {info = 2}  relative error between two consecutive iterates
                      is at most {xtol}.
     
            {info = 3}  conditions for {info = 1} and {info = 2} both hold.
     
            {info = 4}  the cosine of the angle between {fvec} and any
                      column of the Jacobian is at most {gtol} in
                      absolute value.
     
            {info = 5}  number of calls to {fcn} has reached or
                      exceeded {maxfev}.
     
            {info = 6}  {ftol} is too small.  No further reduction in
                      the sum of squares is possible.
     
            {info = 7}  {xtol} is too small.  No further improvement in
                      the approximate solution {x} is possible.
     
            {info = 8}  {gtol} is too small. {fvec} is orthogonal to the
                      columns of the Jacobian to machine precision.
   
      {nfev} is the address of an integer variable, that on output will be
           set to the number of
           calls to fcn.
   
      {fjac} is an output {m} by {n} array.  the upper {n} by {n} submatrix
           of {fjac} contains an upper triangular matrix {r} with
           diagonal elements of nonincreasing magnitude such that
                                  
             {p'*(jac'*jac)*p = r'*r},
   
           where {p} is a permutation matrix and {jac} is the final
           calculated Jacobian. column {j} of {p} is column {ipvt[j]}
           (see below) of the identity matrix.  The lower trapezoidal
           part of {fjac} contains information generated during
           the computation of {r}.
   
      {ldfjac} is a positive integer, not less than {m},
           that specifies the leading dimension of the array {fjac}.
   
      {ipvt} is an integer output array of length {n}, that
           defines a permutation matrix {p} such that {jac*p = q*r},
           where {jac} is the final calculated Jacobian, {q} is
           orthogonal (not stored), and {r} is upper triangular
           with diagonal elements of nonincreasing magnitude.
           Column {j} of {p} is column {ipvt[j]} of the identity matrix.
   
      {qtf} is an output array of length {n} which contains
           the first n elements of the vector {q'*fvec}.
   
      {wa1}, {wa2}, and {wa3} are work arrays of length {n}.
   
      {wa4} is a work array of length {m}.
   
   Other minpack procedures used: {dpmpar,enorm,fdjac2,lmpar,qrfac}.
   
   Argonne National Laboratory, minpack project. March 1980.
   Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More. */

#endif
