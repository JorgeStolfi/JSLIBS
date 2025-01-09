#ifndef sys_stats_H
#define sys_stats_H
#include <stdio.h>

struct process_stats_t{
  //please read man proc
  int pid; //The process ID. 
  char* comm; //The filename of the executable, in parentheses. This is visible whether or not the executable is swapped out. 
  char* state; /*One character from the string "RSDZTW" where R is running,
	  S is sleeping in an interruptible wait,
	  D is waiting in uninterruptible disk sleep,
	  Z is zombie,
	  T is traced or stopped (on a signal),
	  and W is paging. */
  int ppid; //The PID of the parent. 
  int pgrp; //The process group ID of the process.
  int session; //The session ID of the process. 
  int tty_nr; //The tty the process uses. 
  int tpgid; //The process group ID of the process which currently owns the tty that the process is connected to. 
  unsigned long flags; //The kernel flags word of the process. For bit meanings, see the PF_* defines in <linux/sched.h>. Details depend on the kernel version.
  unsigned long minflt; //The number of minor faults the process has made which have not required loading a memory page from disk
  unsigned long cminflt; //The number of minor faults that the process's waited-for children have made.
  unsigned long majflt; //The number of major faults the process has made which have required loading a memory page from disk.
  unsigned long cmajflt; //The number of major faults that the process's waited-for children have made.
  unsigned long utime; //The number of jiffies that this process has been scheduled in user mode, measured  in  clock ticks (divide by sysconf(_SC_CLK_TCK)
  unsigned long stime; //The number of jiffies that this process has been scheduled in kernel mode,, measured  in  clock ticks (divide by sysconf(_SC_CLK_TCK)
  long int cutime; //The number of jiffies that this process's waited-for children have been scheduled in user mode. 
  long int cstime; //The number of jiffies that this process's waited-for children have been scheduled in kernel mode. 
  long int priority; //The standard nice value, plus fifteen. The value is never negative in the kernel. 
  long int nice; //The nice value ranges from 19 (nicest) to -19 (not nice to others). 
  int ZERO; //This value is hard coded to 0 as a placeholder for a removed field. 
  long int itrealvalue; //The time in jiffies before the next SIGALRM is sent to the process due to an interval timer. 
  unsigned long starttime; //The time in jiffies the process started after system boot. 
  unsigned long vsize; //Virtual memory size in bytes. 
  long int rss; //Resident Set Size: number of pages the process has in real memory, minus 3 for administrative purposes. This is just the pages which count towards text, data, or stack space.
  unsigned long rlim; //Current limit in bytes on the rss of the process
  unsigned long startcode; //The address above which program text can run	
  unsigned long endcode; //The address below which program text can run.
  unsigned long startstack; //The address of the start of the stack. 
  unsigned long kstkesp; //The current value of esp (stack pointer), as found in the kernel stack page for the process. 
  unsigned long kstkeip; //The current EIP (instruction pointer).
  unsigned long signal; //The bitmap of pending signals. 
  unsigned long blocked; //The bitmap of blocked signals.
  unsigned long sigignore; //The bitmap of ignored signals.
  unsigned long sigcatch; //The bitmap of caught signals. 
  unsigned long wchan; //This is the "channel" in which the process is waiting. It is the address of a system call,
  unsigned long nswap; //Number of pages swapped (not maintained). 
  unsigned long cnswap;//Cumulative nswap for child processes (not maintained).
  int exit_signal; //Signal to be sent to parent when we die. 
  int processor; //CPU number last executed on. 
  unsigned long rt_priority; //Real-time scheduling priority 
  unsigned long policy; //Scheduling policy
  
};

typedef struct process_stats_t process_stats_t;

process_stats_t get_process_status(void) ;
/* This function returns a stats_t structure with the process informations
   parsed  from /proc/self/stat */

void print_process_status(FILE* arq,process_stats_t st);
/*Prints the content of a {st} into ${arq} in a human readable form*/

double user_cpu_time_usec(void);
/*Returns the time in seconds, including miliseconds*/

#endif