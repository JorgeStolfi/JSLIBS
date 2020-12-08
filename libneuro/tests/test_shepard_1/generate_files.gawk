#! /usr/bin/gawk -f
# Last edited on 2014-01-15 22:58:44 by stolfilocal

BEGIN {
  # Domain interval is {[xlo_xhi]}:
  xlo = 0.0; 
  xhi = 3.0;
  
  # Define the element centers in the interval {[xlo_xhi]}:
  nc = 0;
  split("", ctr); # Center abscissas are {ctr[0..nc-1]}. 
  split("", rad); # Element reference radii are {rad[0..nc-1]}. 
  ctr[nc] = 0.5; rad[nc] = 0.50; nc++;
  ctr[nc] = 1.0; rad[nc] = 0.50; nc++;
  ctr[nc] = 1.8; rad[nc] = 0.30; nc++;
  ctr[nc] = 2.1; rad[nc] = 0.30; nc++;
  # Write the element centers and radii to file {fcenters}:
  fcenters = "out/centers.txt";
  for (ic = 0; ic < nc; ic++) 
    { printf "%10.4f %10.4f\n", ctr[ic], rad[ic] > fcenters; }
  close(fcenters);
  
  # Write weight plotdata to {fweights}, basis plotdta to {fbasis}:
  fweights = "out/weights.txt"
  fbasis = "out/basis.txt"
  nx = 1000;
  for (ix = 0; ix < nx; ix++) 
    { x = xlo + (ix + 0.5)*(xhi - xlo)/nx;
      
      # Evaluate the basic shepard weights {wt[0..nc-1]} at {x}, write to {fweights}:
      split("", wt); 
      printf "%10.6f ", x > fweights;
      for (ic = 0; ic < nc; ic++) 
        { wt[ic] = shepard_weight(ctr[ic], rad[ic], x);
          printf " %10.6f", wt[ic] > fweights;
        }
      printf "\n" > fweights;
      
      # Evaluate the interpolation basis over the interval, write to file {fbasis}:
      split("", bas); 
      eval_basis(x, nc, ctr, wt, bas);
      printf "%10.6f ", x > fbasis;
      for (ic = 0; ic < nc; ic++) { printf " %10.6f", bas[ic] > fbasis; }
      printf "\n" > fbasis;
    }

  close(fbasis);
  close(fweights);
}

function shepard_weight(ct,rd,x,  d,ord,eps,wg,wp)
  { # Shepard radial weight function for element with center {ct} and reference radius {rd}.
    ord = 1;
    eps = 1.0e-10;
    d = (x - ct)/rd;
    if (d < 0) { d = -d; }
    wp = (d == 0.0 ? 1.0/eps : 1.0/(exp(ord*log(d)) + eps));
    wg = exp(-d*d/2);
    return wp*wg;
  }
  
function phi(r,u,z,  d)
  { # The least-squares basis function with index {r} and reference point {u}, 
    # evaluated at {z}.
    if (r == 0)
      { return 1.0; }
    else if (r == 1) 
      { return z-u; }
    else
      { return log(-1); }
  }

function eval_basis(x,nc,ctr,wt,bas,  nb,ic,jc,zi,wi,rb,sb,vr,vs,M,B,N,alpha)
  {
    # Fills {bas[0..nc-1]} with the interpolation basis elements evaluated at {x}.

    nb = 2;
    # Builds an approximation {s[ic](y) = SUM{alpha[rb,ic]*phi[rb](y) : rb \in 0..nb-1}
    # so that {s[ic](ctr[jc])} approximates {ic==jc} in the weighted least squares
    # sense with weight {wt[jc]}. 

    # Builds the least squares system {M*alpha = B} with the specified weights:
    split("", M); # {M[rb,sb]} is {<phi[rb],phi[sb]>}:
    split("", B); # {B[rb,jc]} is {<phi[rb],dat[jc]>} where {dat[jc](z)} is {z == ctr[jc]}:
    for (ic = 0; ic < nc; ic++)
      { # Accumulate the inner product terms for datapoint {ctr[ic]} with weight {wt[ic]}:
        zi = ctr[ic];
        wi = wt[ic];
        for (rb = 0; rb < nb; rb++)
          { vr = phi(rb, x, ctr[ic]);
            for (sb = 0; sb < nb; sb++)
              { vs = phi(sb, x, zi);
                M[rb,sb] += wi*vr*vs; 
              }
            B[rb,ic] = wi*vr;
          }
      }

    # Compute the inverse {N} of matrix {M}:    
    split("", N);
    invert_matrix_2x2(M, N);
    # Multiply {N*B} to get {alpha[0..nb-1,0..nc-1]}:
    split("", alpha);
    multiply_matrices(nb,nb,nc,M,B,alpha);
    # Now evaluate the combination at {x}:
    for (ic = 0; ic < nc; ic++) { bas[ic] = 0; }
    for (rb = 0; rb < nb; rb++)
      { vr = phi(rb, x, x);
        for (ic = 0; ic < nc; ic++)
          { bas[ic] += alpha[rb,ic]*vr; }
      }
  }

function multiply_matrices(m,p,n,A,B,C,  i,j,k,sum)
  {
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { sum = 0;
            for (k = 0; k < p; k++) { sum += A[i,k]*B[k,j]; }
            C[i,j] = sum;
          }
      }
  }

function invert_matrix_2x2(M,N,  det)
  {
    det = M[0,0]*M[1,1] - M[1,0]*M[0,1];
    N[0,0] = M[1,1]/det;
    N[0,1] = -M[0,1]/det;
    N[1,0] = -M[1,0]/det;
    N[1,1] = M[0,0]/det;
  }
