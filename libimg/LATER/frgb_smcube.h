
double frgb_R_from_YUV(double Y, double U, double V);
  /* Computes the smoothed relative saturation {R} of a color {p} given its YUV
    coordinates {Y,U,V}.
    
    The smoothed relative saturation {R} of a color {p} is the linear position of 
    {p} along the segment of colors with the same luminance
    that goes trough {p} and extends from the gray diagonal
    to a point {s} near the boundary of the unit RGB cube. Thus {T} is 0 for a
    gray color, and 1 for any color on the boundary of the RGB
    cube. 
    
    By convention, {T} is zero if {Y} is outside the interval {[0_1]}.
    
    The relative saturation {T} is a continuous function of the RGB
    coordinates, but it is not smooth (C1): it has a kink whenever the
    distal end of the segment crosses an edge of the cube. */

