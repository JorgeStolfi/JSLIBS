#define PROG_NAME "test_sleep"
#define PROG_DESC "test of {jstime.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-17 11:26:39 by stolfi */
/* Created on 2007-01-14 by J. Stolfi, UNICAMP */

#define test_sleep_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include <sys/times.h>
#include <sys/time.h>
/* #include <sys/resource.h> */
/* #include <sys/types.h> */

#include <affirm.h>
#include <jstime.h>
#include <bool.h>

#define MILLION 1000000L
#define BILLION 1000000000L

#define WHICH_MAX 3
  /* Sleeper type {which} is in {0..WHICH_MAX}. */

/* TESTNG FUNCTIONS */

char *get_sleeper_name(uint32_t which);
  /* Returns the name of sleeper procedure number {which}. */
  
double get_sleeper_resolution(uint32_t which);  
  /* Returns the nominal resolution of sleeper procedure {which},
    in microseconds. */
 
double call_sleeper_proc(uint32_t which, double ts);
  /* Calls the sleeper procedure {which} to sleep for {ts} seconds.
    Returns the time that was not slept due to early wakeup, as
    reported by that procedure itself; or {NAN} if the procedure does
    not report it. */

void test_sleeper_proc(uint32_t which, double ts);
  /* Tests the sleeping function selected by {which}, 
    which should be a number in {0..WHICH_MAX},
    with {ts} seconds. */

/* SLEEP FUNCTIONS WITH VARIOUS IMPLEMENTATIONS */

double sleep_nanosleep(clockid_t clk, double ts);
  /* Tries to sleep for {ts} microseconds, using the POSIX
    {clock_nanosleep} function. Returns the difference between {ts} and
    the time actually slept, *as reported by {clock_nanosleep}
    itself*. */

double sleep_usleep(double ts);
  /* Tries to sleep for {ts} seconds, using the BSD
    {usleep} function.  Returns 0. */

double sleep_sleep(double ts);
  /* Tries to sleep for {ts} seconds, using the POSIX {sleep}
    function. If {rus} is not NULL, returns the difference (reported
    by {sleep} itself) between {ts} and the time actually slept. */

double getres_nanosleep(clockid_t clk);
double getres_usleep(void);
double getres_sleep(void);
  /* These functions return the *nominal* granularity of the corresponding 
    sleeper functions, in seconds, as defined in the documentation. */
    
/* MAIN */

int32_t main(int32_t argc, char **argv)
  { fprintf(stderr, "=== testing the various Linux sleeping pills ===\n");
    double maxsleep = 3.000; /* Max duration of sleep attempt (sec). */
    double minsleep = 0.003; /* Min duration of sleep attempt (sec). */ 
    double tsleep = maxsleep;
    while (tsleep >= 0.999999*minsleep)
      { fprintf(stderr, "\n");
        fprintf(stderr, "--- tsleep = %14.12f sec ----------------\n", tsleep); 
        fprintf(stderr, "\n");
        for (int32_t which = 0; which <= WHICH_MAX; which++) 
          { test_sleeper_proc(which, tsleep); }
        tsleep = tsleep/10;
        fprintf(stderr, "\n");
      }
    return 0;
  }

/* ANALYSIS PROCS */

void test_sleeper_proc(uint32_t which, double ts)
  { char *name = get_sleeper_name(which);
    double res = get_sleeper_resolution(which); /* Seconds. */
    fprintf(stderr, "  ... sleeper %d ...\n", which);
    fprintf(stderr, "  %-22s = %-32s\n", "function", name);
    fprintf(stderr, "  %-22s = %28.9f  sec\n", "nominal resolution", res);
    if (ts < res)
      { fprintf(stderr, "  !! not enough resolution -- skipped\n"); }
    else
      { fprintf(stderr, "  %-22s   %28.9f  sec\n", "trying to sleep for", ts);
        double start = real_time_usec()/MILLION;  /* Seconds since epoch. */
        double tleft = call_sleeper_proc(which, ts);
        double stop = real_time_usec()/MILLION;   /* Seconds since epoch. */
        fprintf(stderr, "  %-22s = %28.9f  sec\n", "call time", start);
        if (! isnan(tleft))
          { fprintf(stderr, "  %-22s = %28.9f  sec\n", "reported leftover", tleft); }
        fprintf(stderr, "  %-22s = %28.9f  sec\n", "wakeup time", stop);
        fprintf(stderr, "  %-22s = %28.9f  sec\n", "apparent sleep time", stop - start);
        fprintf(stderr, "  %-22s = %28.9f  sec\n", "apparent sleep error", (stop - start) - ts);
      }
    fprintf(stderr, "\n");
  }

/* SLEEPERS */

double sleep_nanosleep(clockid_t clk, double ts)
  { demand(ts >= 0, "invalid sleep time");
    struct timespec buf, rem;
    uint64_t sec = (uint64_t)floor(ts); /* Integer number of seconds. */
    uint64_t nsec = (uint64_t)floor((ts - (double)sec)*1e9 + 0.5); /* Fraction of second in nanoseconds. */
    /* We should have {0 <= n <= BILLION}, but, just to be sure: */
    while (nsec >= BILLION) { nsec -= BILLION; sec++; }

    assert(sizeof(buf.tv_nsec) == 8);
    buf.tv_nsec = (int64_t)nsec;

    /* fprintf(stderr, "    sizeof(time_t) = %lu\n", sizeof(time_t)); */
    if (sizeof(time_t) == 4)
      { /* This is what the manpage of {time_t} says. */
        assert(sizeof(__time_t) == 4);
        assert(sizeof(buf.tv_sec) == 4);
        demand(sec <= UINT32_MAX, "argument too big");
        buf.tv_sec = (uint32_t)sec;
        fprintf(stderr, "    buf = { %ld, %ld }\n", (int64_t)buf.tv_sec, buf.tv_nsec);
      }
    else if (sizeof(time_t) == 8)
      { /* This is what it is on my Mate/Ubuntu laptop. */
        fprintf(stderr, "    !! {time_t} has non-standard size (64 bits)\n");
        assert(sizeof(__time_t) == 8);
        assert(sizeof(buf.tv_sec) == 8);
        buf.tv_sec = (int64_t)sec;
        fprintf(stderr, "    buf = { %ld, %ld }\n", buf.tv_sec, buf.tv_nsec);
      }
    else
      { demand(FALSE, "bizarre {time_t} size"); }
    rem.tv_sec = rem.tv_nsec = 0; /* Just in case. */
    
    fprintf(stderr, "  %-22s =", "actual request");
    if (sec == 0)
      { fprintf(stderr, " %18s %9ld nsec\n", "", buf.tv_nsec); }
    else
      { fprintf(stderr, " %18ld %09ld nsec\n", buf.tv_sec, buf.tv_nsec); }
    int64_t err = clock_nanosleep(clk, TIMER_ABSTIME, &buf, &rem);
    if (err) 
      { perror("**clock_nanosleep() failed"); 
        fprintf(stderr, "**result = %ld\n", err); 
        exit(1);
      }
    return ((double)rem.tv_sec) + (((double)rem.tv_nsec)/BILLION);
  }

double getres_nanosleep(clockid_t clk)
  { /* The nominal resolution is that of {clock_gettime}: */
    struct timespec buf;
    int32_t err = clock_getres(clk, &buf);
    if (err != 0) 
      { perror("**clock_getres() failed"); 
        fprintf(stderr, "**result = %d\n", err); 
        exit(1);
      }
    return ((double)buf.tv_sec) + (((double)buf.tv_nsec)/BILLION);
  }

double sleep_sleep(double ts)
  { demand(ts >= 0, "invalid sleep time");
    uint64_t sec = (uint64_t)floor(ts + 0.5); /* Integer number of seconds. */
    fprintf(stderr, "  %-22s %20lu            sec\n", "actual request", sec);
    demand(sec <= UINT32_MAX, "argument too big");
    uint32_t rem = sleep((uint32_t)sec);
    return (double)rem;
  }

double getres_sleep(void)
  { /* The nominal resolution is 1 sec: */
    return 1.0;
  }
  
double sleep_usleep(double ts)
  { demand(ts >= 0, "invalid sleep time");
    uint64_t usec = (uint64_t)floor(ts*MILLION + 0.5); /* Integer number of usec. */
    fprintf(stderr, "  %-22s =", "actual request");
    if (usec < MILLION)
      { fprintf(stderr, " %18s %6lu    usec\n", "", usec); }
    else
      { fprintf(stderr, " %18lu %06lu    usec\n", usec/MILLION, usec % MILLION); }
    if (ts >= 1.0) { fprintf(stderr, "  !! warning -- exceeds the official limit (1 sec)\n"); }
    assert(sizeof(useconds_t) == 4);
    demand(usec <= UINT32_MAX, "argument too big");
    int32_t res = usleep((uint32_t)usec);
    if (res < 0)
      { perror("**usleep() failed"); 
        fprintf(stderr, "**result = %d\n", res); 
        exit(1);
      }
    return NAN;
  }

double getres_usleep(void)
  { /* The nominal resolution is 1 microsecond: */
    return 1.0/MILLION;
  }

/* GENERIC SLEEPER FUNCTIONS */

char *get_sleeper_name(uint32_t which)
  { switch(which)
      { case 0: return "nanosleep with CLOCK_MONOTONIC"; break;
        case 1: return "nanosleep with CLOCK_REALTIME"; break;
        case 2: return "POSIX {sleep} function"; break;
        case 3: return "BSD/SUSv2 {usleep} function"; break;
        default: demand(FALSE, "invalid {which} parameter"); return NULL;
      }
  }
 
double get_sleeper_resolution(uint32_t which) 
  { switch(which)
      { case  0: return getres_nanosleep(CLOCK_MONOTONIC);
        case  1: return getres_nanosleep(CLOCK_REALTIME);
        case  2: return getres_sleep();
        case  3: return getres_usleep();
        default: demand(FALSE, "invalid {which} parameter"); return 0;
      }
  }
 
double call_sleeper_proc(uint32_t which, double ts)
  { switch(which)
      { case  0: return sleep_nanosleep(CLOCK_MONOTONIC, ts);
        case  1: return sleep_nanosleep(CLOCK_REALTIME, ts);
        case  2: return sleep_sleep(ts);
        case  3: return sleep_usleep(ts);
        default: demand(FALSE, "invalid {which} parameter"); return 0;
      }
  }
