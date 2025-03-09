# Last edited on 2025-03-03 14:36:18 by stolfi

Data created by Gustavo Lelis for his thesis.

File names are "{yyyy}-{mm}-{dd}-{seq}-{scene}-{method}-{nx:04d}-{ny:04d}-{dtype}.fni"
where {seq} is a two-digit sequence number for same date, 
{dtype} is "Z" for heights, "N" for normals, "G" for gradient, and
{method} is

  "MF" = depth estimation by multifocus stereo.
  "RP" = normal/gradient est by "Robust PCA".
  "WA" = normal/gradient est by Least Squares Woodham.
  "W3" = normal/gradient est by Woodham with 3 brightest images per prixel.
  "SY" = synthetic.

The {scene} may be, among others,

  "manyob" = several balls, disks, pyramids, cones above flat floor plane.