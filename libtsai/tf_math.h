#ifndef tf_math_H
#define tf_math_H

/* Miscellaneous math functions. */

#define SIGNBIT(a) (((a) > 0) ? 0 : 1)

#define SQR(a) ((a) * (a))

#define SQRT3 1.732050807568877293527446341505872366943
#define SQRT(x) sqrt(fabs(x))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CUB(a) ((a)*(a)*(a))


#endif
