/* Timing tools */

#ifndef timefunc_H
#define timefunc_H

#include <stdint.h>

typedef int64_t tfn_func_t(int64_t effort);
  /* The type of a function argument to {time_func}. It should execute
    some other function {f} a certain number of times, proportional to
    {effort}. It should return the actual number of calls to {f}. */

void time_func(char *title, tfn_func_t ex_func, tfn_func_t ex_null, int64_t effort);
  /* Estimates the cpu time per call of a function {f}. Assumes
    that {ex_func(effort)} calls {f} acertain number of times, 
    proportional to [effort}; and {ex_null(effort)} 
    wastes the same loop control overhead, but without 
    actually calling {f}. */

#endif

