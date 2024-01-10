/* See timefunc.h */
/* Last edited on 2011-12-21 20:45:29 by stolfi */

#include <stdio.h>
#include <assert.h>

#include <jstime.h>
#include <jswsize.h>
#include <timefunc.h>

void time_func(char *title, tfn_func_t ex_func, tfn_func_t ex_null, int64_t effort)
  {   
    /* Time an empty loop, to discount... */
    double start = user_cpu_time_usec();
    int64_t ntimes_null = ex_null(effort);
    double stop = user_cpu_time_usec();
    double tare = stop-start;

    /* Now time the function itself: */
    start = user_cpu_time_usec();
    int64_t ntimes_func = ex_func(effort);
    stop = user_cpu_time_usec();
    
    /* They should execute the function the same number of times: */
    assert(ntimes_null == ntimes_func);
    
    /* Compute the times: */
    fprintf(stderr, "times for %s:\n", title);
    fprintf(stderr, ("  total iterations = %" int64_d_fmt "\n"), ntimes_func);
    fprintf(stderr, "  start clock (s) =  %.3f\n", start/1000000.0);
    fprintf(stderr, "  stop clock (s) =   %.3f\n", stop/1000000.0);
    fprintf(stderr, "  total time (s) =   %.3f\n", (stop-start)/1000000.0);
    double tare_pct = 100*tare/(stop-start+0.00001);
    fprintf(stderr, "  overhead (s) =     %.3f", tare/1000000.0);
    fprintf(stderr, " (%.1f%%)\n", tare_pct);
    double avg_time = (stop-start-tare)/((double)ntimes_func);
    fprintf(stderr, "  average ns/call =  %7.1f\n", 1000*avg_time);
  }
 
