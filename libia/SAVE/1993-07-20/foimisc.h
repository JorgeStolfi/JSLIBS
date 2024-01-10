/* miscellaneous definitions */

#ifndef FOIMISC_H
#define FOIMISC_H

#include <stdio.h>

typedef void * MemP;  /* Pointer to memory area */
typedef int MemSize;  /* Size in bytes of memory area */

void error (char *msg);
  /* Prints string $*msg$ to $stderr$ and stops */

void assert(int test, char *msg);
  /* If test is false, prints $*msg$ and stops. */

char *txtcat (char *a, char *b);
  /* Returns a new string that is the concatenation of $*a$ and $*b$. */
  
char *today(void);
  /* Returns today's date in the format "yy-mm-dd hh:mm:ss" */

#endif

