#! /usr/bin/gawk -f 
# Last edited on 2021-08-24 23:55:28 by stolfi

# Creates a FNI image suitable for the value field input of {test_neuromat_image}. 

BEGIN { 
  nx = 280;
  ny = 320; 
  nc = 1; 
  
  rad = 20
  
  printf "begin float_image_t (format of 2006-03-25)\n"; 
  printf "NC = %d\n", nc;
  printf "NX = %d\n", nx;
  printf "NY = %d\n", ny;

  for (iy = 0; iy < ny; iy++) { 
    for (ix = 0; ix < nx; ix++) { 
      printf "%5d %5d", ix, iy;
      for (ic = 0; ic < nc; ic++) { 
        ctrx = nx/2;
        ctry = ny/2;
        dx = ix + 0.5 - ctrx;
        dy = iy + 0.5 - ctry;
        d = sqrt(dx*dx + dy*dy);
        if (d == 0) 
          { v = +1; }
        else
          { z = 6.5*d/rad;
            v = sin(z)/z; 
          }
        printf " %+11.7f", v;
      }
      printf "\n";
    }
    printf "\n";
  }
  printf "end float_image_t\n";
}
  
