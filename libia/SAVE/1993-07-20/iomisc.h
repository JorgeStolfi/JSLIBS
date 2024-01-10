/* Prototype declarations that are missing in stdio.h, time.h: */

#ifndef IOMISC_H
#define IOMISC_H

#include <stdio.h>
#include <sys/time.h>
#include <time.h>

int printf(char *, ...);

int fputc(char, FILE*);
int fputs(char *, FILE*);
int fprintf(FILE *, char *, ...);
int fclose(FILE *);
int fflush(FILE *);
int _flsbuf(unsigned char, FILE*);

time_t time (time_t *);
struct tm *localtime(time_t *);

void srandom(int);
int random(void);
float frandom(void);

#endif
