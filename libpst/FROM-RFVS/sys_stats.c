#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys_stats.h>



process_stats_t get_process_status( void){
  FILE* arq;
  char* name = NULL;
  //int pid = getpid();
  name = "/proc/self/stat" ;
  arq = fopen(name,"r");
  process_stats_t st;
  st.comm = (char*)malloc(sizeof(char)*400);
  st.state = (char*)malloc(sizeof(char)*40);
  fscanf(arq,"%d %s %s %d %d %d %d %d",&st.pid,st.comm,st.state,&st.ppid,&st.pgrp,&st.session,&st.tty_nr,&st.tpgid);
  fscanf(arq,"%lu %lu %lu %lu %lu %lu %lu",&st.flags,&st.minflt,&st.cminflt,&st.majflt,&st.cmajflt,&st.utime,&st.stime);
  fscanf(arq,"%ld %ld %ld %ld %d",&st.cutime,&st.cstime,&st.priority,&st.nice,&st.ZERO);
  fscanf(arq,"%ld %lu %lu %ld %lu",&st.itrealvalue,&st.starttime,&st.vsize,&st.rss,&st.rlim); 
  fscanf(arq,"%lu %lu %lu %lu %lu %lu",&st.startcode,&st.endcode,&st.startstack,&st.kstkesp,&st.kstkeip,&st.signal);
  fscanf(arq,"%lu %lu %lu %lu %lu %lu",&st.blocked,&st.sigignore,&st.sigcatch,&st.wchan,&st.nswap,&st.cnswap);
  fscanf(arq,"%d %d %lu %lu",&st.exit_signal,&st.processor,&st.rt_priority,&st.policy);
  fclose(arq);
  return st;
  
}


void print_process_status(FILE* arq,process_stats_t st){
  fprintf(arq,"PID: %d\ncomm: %s\nstate: %s\nPPDID: %d\nPGRP: %d\nsession: %d\nTTY_nr: %d\nTPGID: %d\n",st.pid,st.comm,st.state,st.ppid,st.pgrp,st.session,st.tty_nr,st.tpgid);
  fprintf(arq,"Flags: %lu\nMinflt: %lu\nCminflt: %lu\nMajflt: %lu\nCmajflt: %lu\nUtime: %lu\nStime: %lu\n",st.flags,st.minflt,st.cminflt,st.majflt,st.cmajflt,st.utime,st.stime);
  fprintf(arq,"Cutime %ld\nCstime: %ld\nPriority: %ld\nNice %ld\nZero: %d\n",st.cutime,st.cstime,st.priority,st.nice,st.ZERO);
  fprintf(arq,"Itrealvalue: %ld\nStartTime: %lu\nVsize: %lu\nRss: %ld\nRlim %lu\n",st.itrealvalue,st.starttime,st.vsize,st.rss,st.rlim); 
  fprintf(arq,"Startcode %lu\nEndcode %lu\nStartStack %lu\nKstkesp: %lu\nKstkeip: %lu\nSignal: %lu\n",st.startcode,st.endcode,st.startstack,st.kstkesp,st.kstkeip,st.signal);
  fprintf(arq,"Blocked: %lu\nSigIgnore: %lu\nSigcatch: %lu\nWchan %lu\nNswap: %lu\nCnswap %lu\n",st.blocked,st.sigignore,st.sigcatch,st.wchan,st.nswap,st.cnswap);
  fprintf(arq,"Exit: %d\nProcessor: %d\nRt_Priority: %lu\nPolicy: %lu\n",st.exit_signal,st.processor,st.rt_priority,st.policy);
}


// double user_cpu_time_usec(void){
// 	struct tms buf;
// 	(void)times(&buf);
// 	return(((double) buf.tms_utime)/((double)sysconf(_SC_CLK_TCK)));
// }
