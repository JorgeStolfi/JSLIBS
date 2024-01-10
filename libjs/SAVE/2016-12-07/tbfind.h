/* Find item in a sorted table by linear interpolation + binary search. */
/* Last edited on 2008-01-16 20:37:50 by stolfi */ 

#ifndef tbfind_H
#define tbfind_H

int tb_find(double f(int i), int iMin, int iMax);
  /* Assumes that {f()} is a function defined over {iMin..iMax}.
    Extends it implicitly with {f(iMin-1) = -1}, {f(iMax+1) = +1}.
    Returns an index {i} in {iMin..iMax+1} such that {f(i-1) < 0} and
    {f(i) >= 0}.
    
    As special cases, if {iMin > iMax} or {f(iMin) >= 0}, returns
    {iMin}, returns {iMin}; else, if {f(iMax) < 0}, returns {iMax+1}.
    In all other cases, returns an index in {iMin+1..iMax} Uses a
    combination of bisection and linear interpolation.
    
    The procedure works with any function {f()}, but is most useful if
    if {f()} is non-decreasing. In that case, if there is some {i} in
    {iMin..iMax} such as {f(i) == 0}, the procedure will return that {i}. */

#endif
