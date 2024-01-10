#ifndef ioprotos_H
#define ioprotos_H

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

/* Prototype declarations that may be missing in stdio.h, time.h: */

int printf(const char *, ...);
int scanf(const char *, ...);

int fputc(int, FILE*);
int fgetc(FILE *);
int fputs(const char *, FILE*);
int fprintf(FILE *, const char *, ...);
int fscanf(FILE *, const char *, ...);
int fclose(FILE *);
int fflush(FILE *);
int _filbuf(FILE *);
long _sysconf(int);

#if (! defined(OSF1V4))
int _flsbuf(unsigned int, FILE*);
#endif

time_t time (time_t *);
struct tm *localtime(const time_t *);

#endif
