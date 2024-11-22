/* See jstime.h */
/* Last edited on 2024-11-16 09:16:14 by stolfi */

#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <affirm.h>

#include <jstime.h>

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
    { int32_t n = (int32_t)strlen(buf);
      char *res = talloc(n+1, char);
      strcpy(res, buf);
      return res;
    }
#undef TODAY_BUFSIZE
  }

/* On systems that have {clock_gettime}, we can use it to get the 
  real time and the user CPU time.  For the pro-rated system CPU time,
  it seems that we must use the {times} function. */

#if (_POSIX_TIMERS)
double real_time_usec(void)
  { struct timespec buf;
    clock_gettime(CLOCK_REALTIME, &buf);
    return (((double)buf.tv_sec)*1000000) + (((double)buf.tv_nsec)/1000);
  }

double user_cpu_time_usec(void) 
  { struct timespec buf;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &buf);
    return (((double)buf.tv_sec)*1000000) + (((double)buf.tv_nsec)/1000);
  }

double system_cpu_time_usec(void) 
  { struct tms buf;
    clock_t t;
    times(&buf);
    t = buf.tms_stime;
    return(1000000.0 * ((double) t)/((double)sysconf(_SC_CLK_TCK)));
  }
#else
/* On systems that do not have {clock_gettime}, we get
  all three times from the {times} function. */

double real_time_usec(void)
  { struct tms buf;
    clock_t etime = times(&buf);
    return(1000000.0 * ((double) etime)/((double)sysconf(_SC_CLK_TCK)));
  }

double user_cpu_time_usec(void) 
  { struct tms buf;
    clock_t etime = times(&buf);
    assert(etime >= 0);
    return(1000000.0 * ((double) buf.tms_utime)/((double)sysconf(_SC_CLK_TCK)));
  }

double system_cpu_time_usec(void) 
  { struct tms buf;
    clock_t etime = times(&buf);
    assert(etime >= 0);
    return(1000000.0 * ((double) buf.tms_stime)/((double)sysconf(_SC_CLK_TCK)));
  }
#endif

