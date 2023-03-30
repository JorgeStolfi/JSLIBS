#define PROG_NAME "test_sleep"
#define PROG_DESC "test of {jstime.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-03-26 11:05:55 by stolfi */
/* Created on 2007-01-14 by J. Stolfi, UNICAMP */

#define test_sleep_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
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

/* SLEEP FUNCTIONS WITH VARIOUS IMPLEMENTATIONS */

double sleep_nanosleep(clockid_t clk, double tus);
  /* Tries to sleep for {tus} microseconds, using the POSIX
    {clock_nanosleep} function. If {rus} is not NULL, returns the
    difference (reported by {clock_nanosleep} itself) between {tus}
    and the time actually slept. */

double sleep_usleep(double tus);
  /* Tries to sleep for {tus} microseconds, using the BSD
    {usleep} function.  Returns 0. */

double sleep_sleep(double tus);
  /* Tries to sleep for {tus} microseconds, using the POSIX {sleep}
    function. If {rus} is not NULL, returns the difference (reported
    by {sleep} itself) between {tus} and the time actually slept. */

double getres_nanosleep(clockid_t clk);
double getres_usleep(void);
double getres_sleep(void);
  /* These functions return the *nominal* granularity of the corresponding 
    sleeper functions, as defined in the documentation. */

/* SLEEP DEBUGGING FUNCTIONS */

#define MAXSLEEPER 3

char *get_sleeper_name(int32_t it);
  /* Returns the name of sleeper procedure number {it}. */
  
double get_sleeper_resolution(int32_t it);  
  /* Returns the nominal resolution of sleeper procedure {it},
    in microseconds. */
  
double call_sleeper_proc(int32_t it, double tus);
  /* Calls the sleeper procedure {it} to sleep for {tus} microseconds.
    Returns the time that was not slept due to early wakeup (as
    reported by that procedure itself). */

void debug_sleeper_proc(int32_t it);
  /* Debugs the sleeping function selected by {it}, 
    which should be a number in {0..MAXSLEEPER}. */
    
/* MAIN */

int32_t main(int32_t argc, char **argv)
  { int32_t it;
    fprintf(stderr, "=== testing the various Linux sleeping pills ===\n");
    for (it = 0; it <= MAXSLEEPER; it++) { debug_sleeper_proc(it); }
    return 0;
  }

/* SLEEPERS */

double sleep_nanosleep(clockid_t clk, double tus)
  { struct timespec buf, rem;
    double ts = tus / MILLION; /* Convert microseconds to seconds. */
    int64_t s = (int32_t)floor(ts); /* Integer number of seconds. */
    int64_t n = (int32_t)floor((ts - (double)s)*1e9 + 0.5); /* Fraction of second in nanoseconds. */
    /* We should have {0 <= n <= BILLION}, but, just to be sure: */
    assert(n >= 0);
    while (n >= BILLION) { n -= BILLION; s++; }
    buf.tv_sec = s;
    buf.tv_nsec = n;
    rem.tv_sec = rem.tv_nsec = 0; /* Just in case. */
    fprintf(stderr, "  %-22s   %17ld.%09ld sec\n", "trying to sleep for", buf.tv_sec, buf.tv_nsec);
    int64_t err = clock_nanosleep(clk, TIMER_ABSTIME, &buf, &rem);
    if (err) 
      { perror("**clock_nanosleep() failed"); 
        fprintf(stderr, "**result = %ld\n", err); 
        exit(1);
      }
    return (((double)rem.tv_sec)*MILLION) + (((double)rem.tv_nsec)/1000);
  }

double getres_nanosleep(clockid_t clk)
  { /* The nominal resolution is that of {clock_gettime}: */
    struct timespec buf;
    int64_t err = clock_getres(clk, &buf);
    if (err) 
      { perror("**clock_getres() failed"); 
        fprintf(stderr, "**result = %ld\n", err); 
        exit(1);
      }
    return (((double)buf.tv_sec)*MILLION) + (((double)buf.tv_nsec)/1000);
  }

double sleep_sleep(double tus)
  { double ts = tus / MILLION; /* Convert microseconds to seconds. */
    uint32_t s = (int32_t)ceil(ts); /* Integer number of seconds. */
    fprintf(stderr, "  %-22s   %17u sec\n", "trying to sleep for", s);
    uint32_t r = sleep(s);
    return ((double)r)*MILLION;
  }

double getres_sleep(void)
  { /* The nominal resolution is 1 sec: */
    return (double)MILLION;
  }
  
double sleep_usleep(double tus)
  { if (tus >= MILLION)
      { fprintf(stderr, "  cannot sleep that much!\n"); }
    else
      { uint32_t us = (int32_t)ceil(tus); /* Integer number of microseconds. */
        fprintf(stderr, "  %-22s   %17u.%06d sec\n", "trying to sleep for", 0, us);
        int32_t res = usleep(us);
        if (res < 0)
          { perror("**usleep() failed"); 
            fprintf(stderr, "**result = %d\n", res); 
            exit(1);
          }
      }
    return 0.0;
  }

double getres_usleep(void)
  { /* The nominal resolution is 1 microsecond: */
    return 1.0;
  }

/* GENERIC SLEEPER FUNCTIONS */

char *get_sleeper_name(int32_t it)
  { switch(it)
      { case 0: return "nanosleep with CLOCK_MONOTONIC"; break;
        case 1: return "nanosleep with CLOCK_REALTIME"; break;
        case 2: return "POSIX {sleep} function"; break;
        case 3: return "BSD/SUSv2 {usleep} function"; break;
        default: demand(FALSE, "invalid {it} parameter"); return NULL;
      }
  }
 
double get_sleeper_resolution(int32_t it) 
  { switch(it)
      { case  0: return getres_nanosleep(CLOCK_MONOTONIC);
        case  1: return getres_nanosleep(CLOCK_REALTIME);
        case  2: return getres_sleep();
        case  3: return getres_usleep();
        default: demand(FALSE, "invalid {it} parameter"); return 0;
      }
  }
 
double call_sleeper_proc(int32_t it, double tus)
  { switch(it)
      { case  0: return sleep_nanosleep(CLOCK_MONOTONIC, tus);
        case  1: return sleep_nanosleep(CLOCK_REALTIME, tus);
        case  2: return sleep_sleep(tus);
        case  3: return sleep_usleep(tus);
        default: demand(FALSE, "invalid {it} parameter"); return 0;
      }
  }

/* ANALYSIS PROCS */

void debug_sleeper_proc(int32_t it)
  { char *name = get_sleeper_name(it);
    double res = get_sleeper_resolution(it);
    fprintf(stderr, "sleeper[%d] = %s\n", it, name);
    fprintf(stderr, "\n");
    fprintf(stderr, "  %-22s = %28.10f us\n", "nominal resolution", res);
    double maxsleep = 3.0e6; /* Max duration of sleep attempt (usec). */
    double minsleep = 3.0e3; /* Min duration of sleep attempt (usec). */ 
    double tsleep = minsleep;
    while (tsleep <= 1.001*maxsleep)
      { fprintf(stderr, "  %-22s   %28.10f us\n", "trying to sleep for", tsleep);
        double start = real_time_usec();
        double tleft = call_sleeper_proc(it, tsleep);
        double stop = real_time_usec();
        fprintf(stderr, "  %-22s = %28.10f us\n", "suspend time", start);
        fprintf(stderr, "  %-22s = %28.10f us\n", "requested sleep", tsleep);
        fprintf(stderr, "  %-22s = %28.10f us\n", "leftover", tleft);
        fprintf(stderr, "  %-22s = %28.10f us\n", "wakeup time", stop);
        fprintf(stderr, "  %-22s = %28.10f us\n", "apparent sleep time", stop - start);
        fprintf(stderr, "  %-22s = %28.10f us\n", "apparent sleep error", (stop - start) - tsleep);
        fprintf(stderr, "\n");
        tsleep *= 10;
      }
  }
