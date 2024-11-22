#define PROG_NAME "test_now"
#define PROG_DESC "test of {gettime_*} and {getres_*} from {jstime.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-18 05:33:32 by stolfi */
/* Created on 2007-01-14 by J. Stolfi, UNICAMP */

#define test_now_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdint.h>
/* #include <sys/types.h> */
/* #include <sys/resource.h> */
#include <sys/time.h>
#include <sys/times.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <bool.h>
#include <jsfile.h>
#include <jstime.h>
#include <affirm.h>

#define MILLION 1000000L
#define BILLION 1000000000L

/* TIMING FUNCTIONS WITH VARIOUS IMPLEMENTATIONS */

double gettime_utctime(void);
  /* Return the value of the unix UTC clock, in
    microseconds, using the POSIX {time} function.  The resolution is
    one second. */

double gettime_nanoclock(clockid_t clk);
  /* Return the value of clock {clk}, in microseconds,
    acccording to the POSIX high-resolution functions
    {clock_gettime} and {clock_getres}. */

double gettime_ansiclock(void);
 /* Return the value of the process's CPU time,
    in microseconds, according to the ANSI {clock} function and the
    {CLOCKS_PER_SECOND} variable. */

double gettime_unixtimes(int32_t mode);
  /* Return a time value, in microseconds, obtained from the UNIX
    {times} function and {sysconf(_SC_CLK_TCK)}. The result depend on
    {mode}:
      
      {mode=0} the CPU time spent by the process in user mode.
      {mode=1} the CPU time spent by the system on behalf of the process.
      {mode=2} the real time elapsed since some fixed epoch.
     
    The epoch for {mode=2} is fixed during the life of the process,
    but may change from process to process. */

double gettime_jstime(int32_t mode);
  /* Return a time value, in microseconds, obtained from 
    the functions in {jstime.h}: 
    
      {mode=0} from {user_cpu_time_usec}.
      {mode=1} from {system_cpu_time_usec}.
      {mode=2} from {real_time_usec}. */

double getres_nanoclock(clockid_t clk);
double getres_ansiclock(void);
double getres_utctime(void);
double getres_unixtimes(int32_t mode);
double getres_jstime(int32_t mode);
  /* These functions return the *nominal* granularity of the corresponding 
    clock functions, as defined in the documentation. */

/* CLOCK DEBUGGING FUNCTIONS */

#define MAXTIMER 11

char *get_timer_name(int32_t it);
char *get_timer_tag(int32_t it);
  /* Returns a descriptive name and a short tag
    for timer procedure number {it}. */
  
double get_timer_value(int32_t it);
  /* Returns the current reading of timer procedure {it},
    converted to microseconds. */
  
double get_timer_resolution(int32_t it);
  /* Returns the nominal resolution of timer procedure {it},
    in microseconds. */

void compute(double *q, double *r, double *s);
  /* The `unit of computation' for the following procedures. It adds 1
    to {q} and twaddles the pair {r,s}. It's safe to call it several
    billion times. */

void test_timer_proc(int32_t it, int64_t nmax);
  /* Debugs the timing function selected by {it}, which should be a
    number in {0..MAXTIMER}. Executes the {compute} operation multiple
    times between clock readings, from a small number up to {nmax}
    times. */
    
void compare_timer_procs(FILE *wr, int32_t nsteps, int64_t nops);
  /* Compares all timing functions by reading them periodically during
    a long computation. Writes to {wr} a table with {MAXTIMER+2}
    colums: the total operation count {n} and {MAXTIMER+1} accumutated
    clock readings (in microseconds). The table will have {nsteps+1}
    rows. The first row should be all zeros, and the procedure
    executes the {compute} operation {nops} times between each row. */

void start_timers(int32_t nt, double beg[]);
void stop_timers(int32_t nt, double beg[]);
  /* These procedures read all timers {0..nt-1} and update {beg[it]}
    appropriately. The sequence
    
      {start_timers(nt,beg); COMPUTE(); stop_timers(nt,beg)}
      
     will put in {beg[it]} the time spent in {COMPUTE()} according to
     timer {it}. The procedures try to discount from {beg[it]} the
     time spent within these procedures themselves (including the
     reading of the clocks). */
  
void plot_timer_table(char *filename);
  /* Plots the timer data table with {gnuplot}, assuming it has been
    written to the file called "{filename}". */
    
/* MAIN */

int32_t main(int32_t argc, char **argv)
  { int32_t it;
  
    bool_t debug_timers = FALSE;
    bool_t compare_timers = TRUE;
    
    if (debug_timers) 
      { fprintf(stderr, "\n=== testing the various Linux clocks ===\n");
        for (it = 0; it <= MAXTIMER; it++)  { test_timer_proc(it, 20*MILLION); }
      }

    if (compare_timers) 
      { fprintf(stderr, "\n=== tabulating the various Linux clocks in parallel ===\n");
        char *filename = "out/test.dat";
        FILE *wr = open_write(filename, TRUE);
        compare_timer_procs(wr, 40, 2*MILLION);
        fclose(wr);
        plot_timer_table(filename);
      }

    fprintf(stderr, "\n=== done ===\n");
    return 0;
  }

void plot_timer_table(char *filename)
  { char *cmd = "gnuplot -persist 1> out/test.got 2> out/test.err";
    /* char *cmd = "cat > out/test.gpl"; */
    FILE *gp = popen(cmd, "w");
    affirm(gp != NULL, "could not start 'gnuplot'");
    fprintf(gp, "set terminal x11\n");
    char *sep = "plot ";
    int32_t it;
    for (it = 0; it <= MAXTIMER; it++)
      { char *tag = get_timer_tag(it);
        fprintf(gp, "%s \\\n", sep);
        fprintf(gp, "  \"%s\" using 1:%-2d title \"%s\" with lp", filename, it+2, tag);
        sep = ",";
      }
    fprintf(gp, "\n");
    fflush(gp);
    sleep(30);
    pclose(gp);
  }

double getres_nanoclock(clockid_t clk)
  { struct timespec buf;
    int64_t err = clock_getres(clk, &buf);
    if (err) 
      { perror("**clock_getres() failed"); 
        fprintf(stderr, "**result = %ld\n", err); 
        exit(1);
      }
    return (((double)buf.tv_sec)*MILLION) + (((double)buf.tv_nsec)/1000);
  }

double gettime_nanoclock(clockid_t clk)
  { struct timespec buf;
    int64_t err = clock_gettime(clk, &buf);
    if (err) 
      { perror("**clock_gettime() failed"); 
        fprintf(stderr, "**result = %ld\n", err); 
        exit(1);
      }
    return (((double)buf.tv_sec)*MILLION) + (((double)buf.tv_nsec)/1000);
  }
  
double getres_utctime(void)
  { return (double)MILLION; }

double gettime_utctime(void)
  { time_t utc = time(NULL);
    if (utc < 0) 
      { perror("**time() failed"); 
        fprintf(stderr, "**result = %ld\n", utc); 
        exit(1);
      }
    return ((double)utc)*MILLION;
  }

double getres_ansiclock(void)
{ return ((double)MILLION)/((double)CLOCKS_PER_SEC); }

double gettime_ansiclock(void)
  { clock_t ticks = clock();
    if (ticks < 0) 
      { perror("**clock() failed"); 
        fprintf(stderr, "**returned = %ld\n", ticks);
        exit(1);
      }
    return ((double)ticks)*MILLION/((double)CLOCKS_PER_SEC);
  }
  
double getres_unixtimes(int32_t mode)
  { return ((double)MILLION)/((double)sysconf(_SC_CLK_TCK)); }

double gettime_unixtimes(int32_t mode)
  { struct tms ustimes;
    clock_t etime = times(&ustimes);
    if (etime < 0) 
      { perror("**times() failed"); 
        fprintf(stderr, "**result = %ld\n", etime); 
        exit(1);
      }
    int64_t CPS = sysconf(_SC_CLK_TCK);
    switch(mode)
      { case 0: 
          return ((double)ustimes.tms_utime)*MILLION/((double)CPS);
        case 1: 
          return ((double)ustimes.tms_stime)*MILLION/((double)CPS);
        case 2: 
          return ((double)etime)*MILLION/((double)CPS);
        default:
          demand(FALSE, "invalid {mode}");
          return 0;
      }
  }

double getres_jstime(int32_t mode)
  { double maxus = 100.0 * 365.0 * 24.0 * 3600.0 * 1.0e6; /* 100 years in microseconds */
    double relres = 1.0e-15; /* Approx relative resolution of a {double}. */
    return maxus*relres; 
  }

double gettime_jstime(int32_t mode)
  { switch(mode)
      {
        case 0: return user_cpu_time_usec(); 
        case 1: return system_cpu_time_usec(); 
        case 2: return real_time_usec(); 
        default: assert(FALSE);
      }
  }
  
/* GENERIC TIMER FUNCTIONS */

char *get_timer_name(int32_t it)
  { switch(it)
      { case  0: return "nanoclock with CLOCK_PROCESS_CPUTIME_ID"; break;
        case  1: return "nanoclock with CLOCK_THREAD_CPUTIME_ID"; break;
        case  2: return "nanoclock with CLOCK_MONOTONIC"; break;
        case  3: return "nanoclock with CLOCK_REALTIME"; break;
        case  4: return "POSIX UTC {time} function"; break;
        case  5: return "ANSI {clock} function"; break;
        case  6: return "Unix {times} function (user)"; break;
        case  7: return "Unix {times} function (system)"; break;
        case  8: return "Unix {times} function (elapsed)"; break;
        case  9: return "Function {user_cpu_time_usec} of {jstime.h}"; break;
        case 10: return "Function {system_cpu_time_usec} of {jstime.h}"; break;
        case 11: return "Function {real_time_usec} of {jstime.h}"; break;
        default: demand(FALSE, "invalid {it} parameter"); return NULL;
      }
  }

char *get_timer_tag(int32_t it)
  { switch(it)
      { case  0: return "nanoc(P)"; break;
        case  1: return "nanoc(T)"; break;
        case  2: return "nanoc(M)"; break;
        case  3: return "nanoc(R)"; break;
        case  4: return "time()"; break;
        case  5: return "clock()"; break;
        case  6: return "times(u)"; break;
        case  7: return "times(s)"; break;
        case  8: return "times(e)"; break;
        case  9: return "jstime(u)"; break;
        case 10: return "jstime(s)"; break;
        case 11: return "jstime(e)"; break;
        default: demand(FALSE, "invalid {it} parameter"); return NULL;
      }
  }

double get_timer_value(int32_t it)
  { switch(it)
      { case  0: return gettime_nanoclock(CLOCK_PROCESS_CPUTIME_ID);
        case  1: return gettime_nanoclock(CLOCK_THREAD_CPUTIME_ID);
        case  2: return gettime_nanoclock(CLOCK_MONOTONIC);
        case  3: return gettime_nanoclock(CLOCK_REALTIME);
        case  4: return gettime_utctime();
        case  5: return gettime_ansiclock();
        case  6: return gettime_unixtimes(0);
        case  7: return gettime_unixtimes(1);
        case  8: return gettime_unixtimes(2);
        case  9: return gettime_jstime(0);
        case 10: return gettime_jstime(1);
        case 11: return gettime_jstime(2);
        default: demand(FALSE, "invalid {it} parameter"); return 0;
      }
  }

double get_timer_resolution(int32_t it)
  { switch(it)
      { case  0: return getres_nanoclock(CLOCK_PROCESS_CPUTIME_ID);
        case  1: return getres_nanoclock(CLOCK_THREAD_CPUTIME_ID);
        case  2: return getres_nanoclock(CLOCK_MONOTONIC);
        case  3: return getres_nanoclock(CLOCK_REALTIME);
        case  4: return getres_utctime();
        case  5: return getres_ansiclock();
        case  6: return getres_unixtimes(0);
        case  7: return getres_unixtimes(1);
        case  8: return getres_unixtimes(2);
        case  9: return getres_jstime(0);
        case 10: return getres_jstime(1);
        case 11: return getres_jstime(2);
        default: demand(FALSE, "invalid {it} parameter"); return 0;
      }
  }

/* ANALYSIS PROCS */

void start_timers(int32_t nt, double beg[])
  { /* We read each timer {it} 4 times, obtaining readings
      {t0,t1,t2,t3}, in such a way that the means {a=(t0+t1)/2} and
      {b=(t2+t3)/2} are coincident for all timers. We then compute
      {D = b + (b-a)/2 = (3b-a)/2 = (3t2+3t3-t0-t1)/4} which 
      should be the correct reading (on the average) when the
      procedure exits.  We set {beg[it]} to {-D} in preparation
      for {stop_timers}. */
    int32_t it;
    for (it = 0; it < nt; it++)    { beg[it] = 0; }
    for (it = 0; it < nt; it++)    { beg[it] += 0.25*get_timer_value(it); }
    for (it = nt-1; it >= 0; it--) { beg[it] += 0.25*get_timer_value(it); }
    for (it = 0; it < nt; it++)    { beg[it] -= 0.75*get_timer_value(it); }
    for (it = nt-1; it >= 0; it--) { beg[it] -= 0.75*get_timer_value(it); }
  }

void stop_timers(int32_t nt, double beg[])
  { /* We do the same trick as in {start_timers}, except that 
    we extrapolate the readings to the time of entering the procedure,
    rather than to its exit.  Namely, we compute 
    {D = a-(b-a)/2 = (3a-b)/2 = (3t0+3t1-t2-t3)}
    and add that to {beg[it]}. */
    int32_t it;
    for (it = 0; it < nt; it++)    { beg[it] += 0.75*get_timer_value(it); }
    for (it = nt-1; it >= 0; it--) { beg[it] += 0.75*get_timer_value(it); }
    for (it = 0; it < nt; it++)    { beg[it] -= 0.25*get_timer_value(it); }
    for (it = nt-1; it >= 0; it--) { beg[it] -= 0.25*get_timer_value(it); }
  }

#define M_SQRT5 2.23606797749978969640


void compute(double *q, double *r, double *s)
  { (*q) += 1.0;
    (*r) = 1/(1 + (*r));
    (*s) += 5/(2*(*r) +1) - M_SQRT5;
  }

void compare_timer_procs(FILE *wr, int32_t nsteps, int64_t nops)
  { int32_t nt = MAXTIMER + 1;
    /* Print timer names: */
    int32_t it;
    for (it = 0; it < nt; it++)
      { char *name = get_timer_name(it);
        double res = get_timer_resolution(it);
        fprintf(wr, "# timer[%02d] res = %28.10f  %s\n", it, res, name);
      }
      
    /* Print header: */
    fprintf(wr, "# %12s", "iterations");
    for (it = 0; it < nt; it++) { fprintf(wr, "    timer[%02d]", it); }
    fprintf(wr, "\n");
    
    fprintf(wr, "# %12s", "------------");
    for (it = 0; it < nt; it++) { fprintf(wr, " %12s", "------------"); }
    fprintf(wr, "\n");
    
    double beg[nt]; /* Clock readings at beginning of iteration. */
    double acc[nt]; /* Accumulated clock readings. */
    for (it = 0; it < nt; it++) { acc[it] = 0; }
        
    /* Read all timers, save them negated in {beg[it]}: */
    start_timers(nt, beg);
    int64_t tops = 0; /* Number of operations performed so far. */
    int32_t st = 0;  /* Number of table steps already done. */
    while(TRUE)
      { /* Read timers and set {beg[it]} to the interval since last {start_timers}: */
        stop_timers(nt, beg);
        /* Accumulate the interval times in {acc}: */
        for (it = 0; it < nt; it++) { acc[it] += beg[it]; }
        /* Print table row: */
        fprintf(wr, "  %12ld", tops);
        for (it = 0; it < nt; it++) { fprintf(wr, " %12.3f", acc[it]); }
        fprintf(wr, "\n");
        /* Are we done? */
        if (st > nsteps) { break; }
        /* Start counting times again: */
        start_timers(nt, beg);
        /* Burn some CPU time: */
        int64_t k;
        double rini = 0.4615; /* A nice number. */
        double r = rini;
        double t = 0.0000;
        double m = 0.0000;
        for (k = 0; k < nops; k++) { compute(&m, &r, &t); }
        /* Make sure that the loop did something: */
        assert(r != rini);
        assert(t != 0);
        assert(m == ((double)nops));
        /* Update operations count and repeat: */
        tops += nops;
        st++;
      }
    fflush(wr);
  }
       
void test_timer_proc(int32_t it, int64_t nmax)
  { char *name = get_timer_name(it);
    double res = get_timer_resolution(it);
    fprintf(stderr, "timer[%02d] res = %28.10f  %s\n", it, res, name);
    fprintf(stderr, "\n");
    int64_t nmin = 100; /* Initial {n}, preferably a power of 10. */
    int64_t pmin = 2;   /* Initial multiplier in each decade (2 or 3). */
    int64_t n = nmin; /* Number of ops already performed. */
    int64_t p = 1; /* Quotient {n/(nmin*10^d)}. */
    while (n <= nmax)
      { fprintf(stderr, "  %-22s = %17ld\n", "iterations", n);
        /* Read the clock: */
        double start = get_timer_value(it);
        /* Do some CPU-intensive computation: */
        int64_t k;
        double rini = 0.4615; /* A nice number. */
        double r = rini;
        double t = 0.0000;
        double m = 0.0000;
        for (k = 0; k < n; k++) { compute(&m, &r, &t); }
        /* Read the clock again: */
        double stop = get_timer_value(it);
        /* Make sure that the loop did something: */
        assert(r != rini);
        assert(t != 0);
        assert(m == ((double)n));
        /* Print results: */
        fprintf(stderr, "  %-22s = %28.10f us\n", "start time", start);
        fprintf(stderr, "  %-22s = %28.10f us\n", "stop time", stop);
        fprintf(stderr, "  %-22s = %28.10f us\n", "time interval", stop - start);
        fprintf(stderr, "  %-22s = %28.10f us\n", "time per iteration", (stop - start)/((double)n));
        fprintf(stderr, "\n");
        /* Go to the next value: */
        n = n/p;
        if (p == 1)
          { /* Start of a new decade: */ p = pmin; }
        else if (p == 2)
          { /* Sequence: 1,2,5,10,20,50,100,200, ... */ p = 5; }
        else if (p == 3)
          { /* Sequence: 1,3,10,30,100,300,... */ n = 10*n; p = 1; }
        else if (p == 5)
          { /* Sequence: 1,2,5,10,20,50,100,200, ... */ n = 10*n; p = 1; }
        else
          { assert(FALSE); }
        n = p*n;
      }
  }
