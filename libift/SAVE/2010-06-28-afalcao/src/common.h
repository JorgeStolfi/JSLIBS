#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>

/* Error messages */

#define MSG1  "Cannot allocate memory space in"
#define MSG2  "Cannot open file in"

/* Common data types to all programs */ 

typedef enum{false,true} bool; 
typedef struct timeval timer;
typedef unsigned char uchar;

/* Common definitions */


#define PI          3.1415926536
#define INTERIOR    0
#define EXTERIOR    1
#define BOTH        2
#define WHITE       0 
#define GRAY        1
#define BLACK       2
#define NIL        -1
#define INCREASING  1
#define DECREASING  0
#define Epsilon     1E-05       

/* Common operations */

#define MAX(x,y) ((x > y)?x:y)
#define MIN(x,y) ((x < y)?x:y)
#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))

char   *AllocCharArray(int n);  /* It allocates 1D array of n characters */
uchar  *AllocUCharArray(int n);  /* It allocates 1D array of n characters */
ushort *AllocUShortArray(int n);  /* It allocates 1D array of n characters */
int    *AllocIntArray(int n);   /* It allocates 1D array of n integers */
float  *AllocFloatArray(int n); /* It allocates 1D array of n floats */
double *AllocDoubleArray(int n);/* It allocates 1D array of n doubles */

void Error(char *msg,char *func); /* It prints error message and exits
                                     the program. */
void Warning(char *msg,char *func); /* It prints warning message and
                                       leaves the routine. */
void Change(int *a, int *b); /* It changes content between a and b */

timer *Tic(); /* It marks the initial time */
timer *Toc(); /* It marks the final time */
float CTime(timer *tic, timer *toc); /* It computes the time difference */

#endif








