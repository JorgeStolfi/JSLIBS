/* Last edited on 2009-03-02 22:07:35 by stolfi */

typedef struct pst_camera_t 
  { hr3_point_t O;  /* The view point of the camera. */
  } pst_camera_t;
  /* Representation of a camera. The viewpoint is given as a quadruple
    {O} of signed homogeneous coordinates {[m,x,y,z]}, where {m}
    is the /weight/ (scale factor). The {z} coordinate is always
    positive, and {m} is always non-negative. If {m} is zero, the
    viewpoint is at infinity in the direction of the vector {(x,y,z)}.
    If {m} is positive, the cartesian coordinates {(O.X,O.Y,O.Z)} of
    the viewpoint are {(x/m,y/m,z/m)}. */

