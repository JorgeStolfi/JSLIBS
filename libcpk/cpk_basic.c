/* See cpk_basic.h */
/* Last edited on 2024-12-31 14:36:00 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
/* #include <sys/resource.h> */
/* #include <sys/types.h> */
#include <sys/times.h>
#include <limits.h>

#include <affirm.h>
#include <jsstring.h>
#include <vec.h>
#include <bool.h>

#include <cpk_basic.h>

char *today(void)
  {
#define TODAY_BUFSIZE 200
    char buf[TODAY_BUFSIZE];
    int32_t rcode;
    time_t today_secs = time(NULL);
    struct tm today;
    today = *localtime(&today_secs);
    rcode = snprintf(buf, TODAY_BUFSIZE,
      "%02d-%02d-%02d %02d:%02d:%02d", 
      today.tm_year % 100, today.tm_mon, today.tm_mday, 
      today.tm_hour, today.tm_min, today.tm_sec
    );
    affirm (rcode >= 0, "snprintf failed");
    { size_t n = strlen(buf);
      char *res = (char *)malloc(n+1);
      strcpy(res, buf);
      return res;
    }
#undef TODAY_BUFSIZE
  }

/* extern long int32_t __sysconf (int32_t); */
double cpk_cpu_time_1(void) 
  { struct tms buf;
    clock_t t;
    times(&buf);
    t = buf.tms_utime;
    return(1000000.0 * ((double) t)/((double) sysconf(_SC_CLK_TCK)));
  }

double cpk_cpu_time_2(void) 
  { struct timespec buf;
    clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &buf);
    return (((double)buf.tv_sec)*1000000) + (((double)buf.tv_nsec)/1000);
  }

int32_t intcmp(int32_t x, int32_t y) 
  { return (x < y ? -1 : ( x > y ? +1 : 0)); }
  
int32_t dblcmp(double x, double y) 
  { return (x < y ? -1 : ( x > y ? +1 : 0)); }

vec_typeimpl(ui2_vec_t,ui2_vec,ui2_t);
