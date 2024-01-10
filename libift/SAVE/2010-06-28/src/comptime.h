#ifndef _CTIME_H_
#define _CTIME_H_

#include "common.h"

timer *Tic(); /* It marks the initial time */
timer *Toc(); /* It marks the final time */
float CTime(timer *tic, timer *toc); /* It computes the time difference */

#endif
