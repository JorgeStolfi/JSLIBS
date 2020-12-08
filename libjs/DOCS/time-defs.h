/*  Last edited on 2007-01-13 17:38:35 by stolfi */

/* A collection of time-related types, constants, and functions extracted from the 
  GNU /usr/include/*.h mess. Beware that some of these symbols will be available
  only in some dialects of C.  */

/* TYPES */

  /* Time value in CPU clock ticks: */
  typedef long int clock_t;

  /* Time value, usually in seconds (or ticks?): */
  typedef long int time_t;

  /* Fraction of a second in microseconds: */
  typedef long int suseconds_t;
  
  /* Time value, in nanoseconds: */
  struct timespec
    { time_t tv_sec;    /* Seconds.  */
      long int tv_nsec; /* Nanoseconds; must be in {0..999999999}.  */
    };

  /* Time value, in microseconds: */
  struct timeval
    { time_t tv_sec;       /* Seconds.  */
      suseconds_t tv_usec; /* Microseconds.  */
    };

  /* POSIX.1b structure for timer start values and intervals: */
  struct itimerspec
    { struct timespec it_interval;
      struct timespec it_value;
    };

  /* Time value, broken down into conventional components: */
  struct tm
    { int tm_sec;           /* Seconds.     [0-60] (1 leap second) */
      int tm_min;           /* Minutes.     [0-59] */
      int tm_hour;          /* Hours.       [0-23] */
      int tm_mday;          /* Day.         [1-31] */
      int tm_mon;           /* Month.       [0-11] */
      int tm_year;          /* Year - 1900.  */
      int tm_wday;          /* Day of week. [0-6] */
      int tm_yday;          /* Days in year.[0-365] */
      int tm_isdst;         /* DST.         [-1/0/1]*/
      long int tm_gmtoff;   /* Seconds east of UTC.  */
      char *tm_zone;        /* Timezone abbreviation.  */
    };
    /* The fields {tm_gmtoff} and {tm_zone} are available only if {__USE_BSD}
      is defined.  In other systems, they are called {__tm_gmtoff} and {__tm_zone}. */

  /* Structure describing CPU time used by a process and its children.  */
  struct tms
    { clock_t tms_utime;    /* User CPU time.  */
      clock_t tms_stime;    /* System CPU time.  */
      clock_t tms_cutime;   /* User CPU time of dead children.  */
      clock_t tms_cstime;   /* System CPU time of dead children.  */
    };
 
/* FUNCTIONS */

  time_t time(time_t *t);
    /*  Return the current time since the Epoch (00:00:00 UTC, January 1,
      1970), measured in seconds. If {t} is not NULL, also puts it in {*t}.
      Returns -1 and sets {errno} if the clock is not available. 

      POSIX.1  defines seconds since the Epoch as a value to be interpreted as the number of
      seconds between a specified time and the Epoch, according to a formula for  conversion
      from  UTC  equivalent  to conversion on the na??ve basis that leap seconds are ignored
      and all years divisible by 4 are leap years.  This value is not the same as the actual
      number  of seconds between the time and the Epoch, because of leap seconds and because
      clocks are not required to be synchronised to a standard reference.  The intention  is
      that  the  interpretation of seconds since the Epoch values be consistent; see POSIX.1
      Annex B 2.2.2 for further rationale.

      Conforming to: SVr4, SVID, POSIX, X/OPEN, 4.3BSD.
      Under 4.3BSD, this call is obsoleted by gettimeofday(2). */

   clock_t clock(void);
    /* Returns the processor time used by the program so far (user
      time + system time).  Returns -1 if
      the clock is not available.
      
      The result is in a special unit; it should be divided by
      {CLOCKS_PER_SECOND} to yield the program time in seconds

      Linux does not include the times of waited-for children in the
      value returned by {clock}. The {times} function, which
      explicitly returns (separate) information about the caller and
      its children, may be preferable.

      Conforms to: ANSI  C. */

 clock_t times(tms *buf);
    /* Stores into {*buf} the CPU time used by this process and all its
      dead children (and their dead children) in BUFFER.
      
      On success, return the elapsed real time from some moment in the
      past, possibly dependent on system boot time. This may overflow
      the range of a {clock_t}. Returns -1 for errors. 
      
      All times are in clock ticks, which are different from those of
      the {clock} function. The number of clock ticks per second can
      be obtained using {sysconf(_SC_CLK_TCK)}. The symbols {HZ} and
      {CLK_TCK} of older systems are now obsolete.

      Conforms to: SVr4, SVID, POSIX, X/OPEN, 4.3BSD. */
  
  double difftime(time_t t1, time_t t0);
    /* Returns the difference {t1 - t0}, in seconds.
    
      Conforming to: SVID 3, 4.3BSD, ISO 9899. */

  time_t mktime(tm *tp);
  time_t timelocal(tm *tp);
    /* These functions are synonymous. They
      convert a broken-down time structure {*tp}, expressed as local time,
      to calendar time representation.  The function ignores the
      current values of {tp->tm_wday} and {tp->tm_yday} and recomputes them from 
      the other {tp} fields.  If any {tp} fields are  outside  their  legal
      interval,  they will be normalized (so that, e.g., 40 October is changed into 9 Novem-
      ber).  Calling {mktime()} also sets the external variable {tzname} with information  about
      the  current  time  zone.   If the contents of {*tp} cannot be represented as
      calendar time (seconds since the epoch), {mktime()} returns -1  and
      does not alter the {tp->tm_wday} and {tp->tm_yday}.

      Conforming to: {mktime} in SVID 3, POSIX, 4.3BSD, ISO 9899; {timelocal} is a 
      nonstandard GNU extension. */

  time_t timegm(tm *tp);
    /*  Equivalent {mktime} but always interprets {*tp} as a
      Coordinated Universal Time (UTC), not local time. Equivalent to
      saving the (possibly NULL) value of the environment variable
      {TZ}, setting it to "", calling {mktime}, and restoring the
      value of {TZ}.

      Conforming to: None (it's a nonstandard GNU extension).  */

  size_t strftime(char *s, size_t maxsize, char *fmt, tm *tp);
    /* Format {*tp} into {*s} according to format {*fmt}.
      See the manpage for the '%'-conversions allowed in {fmt}.
      Write no more than {maxsize} characters and return the number
      of characters written, or 0 if it would exceed {maxsize}. 
      
      Conforming to: ANSI  C, SVID 3, ISO 9899.  Some '%'-conversions
      are standard-dependent. */

  char *strptime(char *s, char *fmt, tm *tp);
    /* The converse of {strftime}: it parses the string {*s} according
      to {*fmt} and stores the corresponding binary time
      information in {*tp}. See the manpage for the '%'-conversions
      allowed in {fmt}. The return value is a pointer to the first
      unparsed character in {*s} (possibly the zero byte of {*s}).
      Returns NULL in case of error. 
      
      Conforming to: XPG4, SUSv2, POSIX 1003.1-2001. */

  size_t strftime_l(char *s, size_t maxsize, char *fmt, tm *tp, locale_t loc);
  char *strptime_l(char *s, char *fmt, tm *tp, locale_t loc);
    /* Similar to {strftime,strptime} but take the information from
      the provided locale and not the global locale.  */

  tm *gmtime(time_t *t);
  tm *gmtime_r(time_t *t, tm *tp);
    /* These functions break down the time {*t} (interpreted as
      seconds since the Epoch) into its conventional components,
      according to the Coordinated Universal Time (UTC) standard and
      returns a pointer to the resulting broken-down time structure.
      The procedures return NULL in case of failure.
      
      The return value of {gmtime} points to a statically allocated
      struct which may be overwritten by other date/time functions, and is not
      thread-safe. The {gmtime_r} stores the data in the
      client-supplied struct {*tp} and returns {tp}.
      
      Conforming to: SVID 3, POSIX, 4.3BSD, ISO 9899. The {localtime_r} version
      is in SUSv2. */

  tm *localtime(time_t *t);
  tm *localtime_r(time_t *t, tm *tp);
    /* These functions break down the time {*t} (interpreted as
      seconds since the Epoch) into its conventional components,
      according to the local time zone and daylight savings. 
      The procedures return NULL in case of failure.

      The function {localtime} also sets the global variables
      {tzname,timezone,daylight} as explained under {tzset}. Its return
      value points to a statically allocated struct which may be
      overwritten by other date/time functions and is not thread-safe.

      The {localtime_r} function stores its result in a user-supplied
      area {*tp}. Apparently (???), this function does not modify
      {tzname,timezone,daylight}; instead, it uses (or sets?) the fields
      {tp->tm_gmtoff,tp->tm_zone}.

      Conforming_to: SVID 3, POSIX, 4.3BSD, ISO 9899. The {localtime_r} version
      is in SUSv2. */

  char *asctime(tm *tp);
  char *asctime_r(tm *tp, char *buf);
    /* These procedures convert the broken-down time {*tp} to a string
      of the form "{Wdy} {Mth} {dd} {hh}:{mm}:{ss} {yyyy}\n" and
      return a pointer to that string. The procedures return NULL in
      case of failure.
     
      The function {asctime} puts the string into a statically
      allocated buffer, which may be overwritten by other date/time
      functions and is not thread-safe. The {asctime_r} function
      stores the result in a user-supplied buffer {buf}, whose length
      must be at least 26.
      
     Conforming to: SVID 3, POSIX, 4.3BSD, ISO 9899. The {asctime_r} version
     is in SUSv2.  */

  char *ctime(time_t *t);
    /* Equivalent to {asctime(localtime(t))}. */
  
  char *ctime_r(time_t *t, char *buf);
    /* Equivalent to {asctime_r(localtime_r(timer, tmp), buf)},
      where {tmp} is an internally allocated temporary {tm} area. */

  char *tzname[2];   /* Current timezone names. */
  long int timezone; /* Seconds west of UTC.  */
  int daylight;      /* If daylight-saving time is ever in use. */
    /* These variables are set by {tzset} below (and some other functions
      that call {tzset} as a side effect). */ 

  void tzset(void);
    /* Sets the variables {tzname,timezone,daylight} from the {TZ}
      environment variable. See the manpage for the valid {TZ} syntax.
      If {TZ} is not defined, a locale-dependent default (read from
      the appropriate system configuration file) is used. If the {TZ}
      environment variable is set but its value is empty or not a
      valid timezone spec, the variables {tzname,timezone,daylight}
      are set for Coordinated Universal Time (UTC).
      
      Conforming to: SVID 3, POSIX, 4.3BSD.  The variable {daylight}
      is obsolete now, but is still required by SUSv2. */

  int dysize(int year);
    /* Return the number of days in {year}. The result is 366 if the
      year is divisible by 4 but not by 100, or is divisible by 400.
      Note that it does not work for years before the Gregorian
      Calendar reform.
      
      Conforming to: SunOS 4.x. */

  unsigned int sleep(unsigned int s);
    /* Suspends the current process until {s} seconds have elapsed, or a signal
      arrives which is not ignored.  Returns 0 if the requested time has elapsed,
      otherwise returns the number of seconds left to sleep.
    
      Conforming to: POSIX.1. */
      
  void usleep(unsigned long u); /* BSD version. */
  int usleep(useconds_t u);     /* SUSv2 version. */
    /* Suspends the current process until {u} microseconds have elapsed, or a signal
      arrives which is not ignored.  
      
      The BSD version returns no value. The SUSv2 version returns
      0 if the sleep was complete, otherwise returns -1 and sets {errno} to {EINTR}
      (interrupted by signal) or {EINVAL} (on some systems, if {u > 1000000}).
    
      Conforming to: BSD, SUSv2. */

  int nanosleep(timespec *req, timespec *rem);
    /* This function suspends the current thread by the amount of time
      specified in {req} or more. In case of success, returns 0.
     
      The function can return earlier if a signal has been delivered
      to the process. In this case, it returns -1, sets {errno} to
      {EINTR}, and saves the remaining time into {*rem} unless {rem}
      is NULL. The value of {*rem} can then be used to call
      {nanosleep} again and complete the specified pause. It also
      returns -1 if {*req} is invalid, setting {errno} to {EFAULT} or
      {EINVAL}.
      
      Compared to {sleep} and {usleep}, {nanosleep} has the advantages
      of not affecting any signals, being POSIX-standard, providing
      higher timing resolution, and making it easier to resume the
      pause.
      
      On the other hand, the implementation of {nanosleep} in Linux
      2.6 is based on the normal kernel timer mechanism, which has a
      resolution of {1/HZ} seconds --- namely, 10 ms on Linux/i386 and
      1 ms on Linux/Alpha. Therefore, {nanosleep} may take {1/HZ}
      seconds longer than specified by {*req}; and {*rem} is usually
      rounded to the next larger multiple of {1/HZ} sec.
      
      Conforming to: POSIX.1b. */

  /* Identifier for system-wide realtime clock.  */
  #define CLOCK_REALTIME             0
  /* Monotonic system-wide clock.  */
  #define CLOCK_MONOTONIC            1
  /* High-resolution timer from the CPU.  */
  #define CLOCK_PROCESS_CPUTIME_ID   2
  /* Thread-specific CPU-time clock.  */
  #define CLOCK_THREAD_CPUTIME_ID    3

  /* Flag to indicate time is absolute.  */
  #define TIMER_ABSTIME  1

  int clock_gettime(clockid_t clock_id, timespec *tp);
  int clock_getres(clockid_t clock_id, timespec *res);
    /* The function {clock_gettime} saves in {*res} the current value
      of clock {clock_id}, whereas {clock_getres} saves that clock's
      nominal resolution. 
      
      Both functions return 0 on success; otherwise, they return
      {-EINVAL} for an invalid {clock_id}, and {-EFAULT} for an
      invalid {tp}.
      
      Conforming to: POSIX.1b.  See under {nanosleep}
      above the caveats for the Linux implementation. */

  int clock_settime(clockid_t clock_id, timespec *tp);
    /* Sets the clock {clock_id} to the value {*tp}.
      
      Conforming to: POSIX.1b.  See under {nanosleep}
      above the caveats for the Linux implementation. */

  int clock_nanosleep(clockid_t clock_id, int flags, timespec *req, timespec *rem);
    /* Like {nanosleep}, but uses the clock {clock_id}.
    
      Returns 0 if the time specified by {*req} has elapsed; otherwise,
      {-EINVAL} for an invalid {clock_id}, and {-EFAULT} for an
      invalid {tp}.
      
      Conforming to: POSIX.1b.  See under {nanosleep}
      above the caveats for the Linux implementation. */

  int clock_getcpuclockid(pid_t pid, clockid_t *clock_id);
    /* Returns in {*clock_id} the ID for CPU-time clock of the 
      process with ID {pid} (or of the calling process, if {pid} is
      zero).
    
      Returns 0 on success; otherwise, returns {EPERM} if the caller
      does not have permission to read the clok of process {pid}; and
      {ESRCH} if there is no such process. */

  int timer_create(clockid_t clock_id, sigevent *evp, timer_t *timer_id);
    /* Creates a new per-process timer based on clock {clock_id}, and
      returns its ID in {*timer_id}. The {*evp} structure defines the
      asynchronous notification that occurs when the timer expires.
      
      Returns 0 on success; otherwise, returns {-EINVAL} (invalid
      {clock_id}), {-EAGAIN} (failed to create the timer), {-EFAULT}
      (bad {*evp}).
      
      Conforming to: POSIX.1b. See under {nanosleep} above the caveats
      for the Linux implementation. */

  int timer_delete(timer_t timer_id);
    /* Deletes timer {timer_id}. Returns 0 on success; otherwise, returns {-EINVAL} (invalid
      {timer_id}). */

  int timer_settime(timer_t timer_id, int flags, itimerspec *v, itimerspec *ov);
    /* If {v->it_value} is not zero, sets timer {timer_id} to expire in
      {v->tv_value}, and starts it running (if not running already).
      If {v->it_value} is zero, the timer is stopped. If {ov} is not
      NULL, also saves the old value in {*ov} (with {ov->it_value}
      zero if the timer was disarmed).
      
      Returns 0 on success; otherwise, returns {-EINVAL} (invalid
      {timer_id}) or {-EFAULT} (invalid {*v} or {*ov}). */

  int timer_gettime(timer_t timer_id, itimerspec *v);
    /* Gets the current value of timer {timer_id} and stores it in {*v}.
      
      Returns 0 on success; otherwise, returns {-EINVAL} (invalid
      {timer_id}) or {-EFAULT} (invalid {*v}).  */

  int timer_getoverrun(timer_t timer_id);
    /* Get the expiration overrun for timer {timer_id}. */


  int getdate_err;
    /* This variable is set by {getdate} to one of the following
      values to indicate an error.
    
       1  The DATEMSK environment variable is null or undefined.
       2  The template file cannot be opened for reading.
       3  Failed to get file status information.
       4  The template file is not a regular file.
       5  An error is encountered while reading the template file.
       6  Memory allication failed (not enough memory available).
       7  There is no line in the template that matches the input.
       8  Invalid input specification Example: February 31 or a time is
          specified that can not be represented in a {time_t}. */

  tm *getdate(char *s);
  int getdate_r (char *s, tm *tp);
    /* These functions parse the string {*s} as a date specification,
      and return it as a broken-down time structure. Unlike
      {strptime}, which takes a format argument, these functions use the
      formats found in the file of which the full path name is given
      in the environment variable {DATEMSK}. The first line in the
      file that matches the given input string is used for the
      conversion.
      
      The function {getdate} is not thread-safe, since it sets the
      variable {getdate_err} in case of error, and stores the result
      in a statically allocated buffer, which may be overwritten by
      other date/time functions. The {getdate_r} function stores the
      data in the client-supplied struct {*tp}; it returns 0 on
      success, and one of the above error codes otherwise. */

/* ---------------------------------------------------------------------- */
/* TIMERS */

  #define ITIMER_REAL  ???
    /* This timer decrements in real time, and delivers {SIGALRM} upon expiration. */

  #define ITIMER_VIRTUAL ???
    /* This timer decrements  only  when the process is executing, and delivers
      {SIGVTALRM} upon expiration. */

  #define ITIMER_PROF ???
    /* This timer decrements both when the
      process executes and when the system is executing on behalf of
      the process. Coupled with {ITIMER_VIRTUAL}, this timer is
      usually used to profile the time spent by the application in
      user and kernel space. Delivers {SIGPROF} upon expiration. */

  int getitimer(int which, struct itimerval *value);
  int setitimer(int which, const struct itimerval *value, struct itimerval *ovalue);
