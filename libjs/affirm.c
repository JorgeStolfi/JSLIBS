/* See affirm.h */
/* Last edited on 2022-10-18 21:09:50 by stolfi */

#include <affirm.h>
#include <stdio.h>
#include <stdlib.h>

void programerror (const char *msg, const char *file, unsigned int line, const char* proc)
  { fprintf (stderr, "%s:%u: ** (%s) %s\n", file, line, proc, msg);
    exit(1);
  }

void *checknotnull(void *p, const char *msg, const char *file, unsigned int line, const char* proc)
  { if (p == NULL) { programerror(msg, file, line, proc); }
    return (void *)p;
  }
