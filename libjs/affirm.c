/* See affirm.h */
/* Last edited on 2019-04-25 17:58:13 by jstolfi */

#include <affirm.h>
#include <stdio.h>
#include <stdlib.h>

void programerror (const char *msg, const char *file, unsigned int line, const char* proc)
  { fprintf (stderr, "%s:%u: ** (%s) %s\n", file, line, proc, msg);
    exit(1);
  }

void *checknotnull(const void *p, const char *msg, const char *file, unsigned int line, const char* proc)
  { if (p == NULL) { programerror(msg, file, line, proc); }
    return (void *)p;
  }
