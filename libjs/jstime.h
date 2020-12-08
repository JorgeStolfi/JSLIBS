#ifndef jstime_H
#define jstime_H

/* J. Stolfi's miscellaneous time utilities. */
/* Last edited on 2009-10-27 18:36:54 by stolfi */

char *today(void);
  /* Returns a newly allocated string containing 
    today's date in the format "yy-mm-dd hh:mm:ss" */

/* TIMING FUNCTIONS

  These functions return a time value measured in
  microseconds and packaged as a {double}. 
  Beware that the granularity of the clock 
  may be greater than (or less than) 1 usec. */

double real_time_usec(void);
  /* Returns the actual time in microseconds since
    some moment in the past.  The actual value is 
    not meaningful, only the differences.*/

double user_cpu_time_usec(void);
double system_cpu_time_usec(void);
  /* Returns the accumulated CPU time of the process, in microseconds:
    in user-mode computation and in system calls on behalf of the
    process, respectively. */

#endif
