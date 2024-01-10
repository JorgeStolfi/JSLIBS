/* miscellaneous functions */

#include "foimisc.h"
#include "iomisc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

void error (char *msg)
  { fprintf (stderr, "*** %s\n", msg);
    exit(1);
  }

void assert(int test, char *msg)
  {
    if (! test)
      { fprintf(stderr, "\n*** assertion is false: %s ***\n", msg);
        exit(1);
      }
  }

char *txtcat (char *a, char *b)
  {
    char *r = malloc(strlen(a)+strlen(b)+1);
    if (r)
      { strcpy(r, a);
        strcat(r, b);
        return(r);
      }
    else
      { error ("textcat: memory exhausted"); return(NULL); }
  }

char *today(void)
  {
    time_t today_secs = time(NULL);
    struct tm today = *localtime(&today_secs);
    char *buf = (char *) malloc(20);
    sprintf(buf, 
      "%02d-%02d-%02d %02d:%02d:%02d\n", 
      today.tm_year % 100, today.tm_mon, today.tm_mday, 
      today.tm_hour, today.tm_min, today.tm_sec
    );
    return(buf);
  }
