#! /usr/bin/gawk -f
# Last edited on 2014-07-29 01:48:48 by stolfilocal

BEGIN \
  { abort = -1;
    usage = ( "subsampling-weights \\\n" \
      "  -v balanced=BOOL \\\n" \
      "  -v variance=NUM \\\n" \
      "  -v radius=NUM \\\n" \
      "  > OUTFILE.wts" \
    );

    # Computes filter weights useful for subsampling
    # at 1/2 the original sampling frequency.

    if (balanced == "") { arg_error("must define \"balanced\""); }
    if (variance == "") { arg_error("must define \"variance\""); }
    if (radius == "") { arg_error("must define \"radius\""); }
    
    n = 2*radius+1;

    split("", owt); # Initial (Gaussian) weights.
    split("", swt); # Balanced, normalized, scaled and quantized weights.
    
    # Adjust the initial deviation {dev} and scale factor {fac}.
    maxfac = 9999.0; # Max scaling factor.
    Nvar = 200;      # Number of {var} values to try on each side of {variance}.
    Nfac = 200;      # Number of {fac} values to try.
    Dfac = 24;       # Denominator of {fac}.
    bvar = -9999;    # True variance that got closer to {variance}.
    bxvr = -1;       # Best value of {xvr}.
    bfac = -1;       # Best value of {fac}. 
    n_accept = 0;    # Count valid candidates.
    n_reject = 0;    # Count invalid candidates.
    for (tvar = 0; tvar <= 2*Nvar; tvar++)
      { # Compute {kvar = 0,-1,+1,-2,+2,-3,+3,...-Nvar,+Nvar} 
        kvar = int((tvar+1)/2)*(1 - 2*(tvar%2));
        # Pick a nominal variance {xvr} near {variance}
        xvr = variance*exp(kvar*log(2.0)/Nvar)
        # printf "xvr = %s\n", xvr > "/dev/stderr";
        for (tfac = 0; tfac < Nfac; tfac++)
          { fac = maxfac - tfac/Dfac;
            # Compute table with parameters {var,fac}:
            init_weights(xvr,n,owt);
            normalize_weights(n,fac,owt,swt,balanced);
            if ((! balanced) || (weight_unbalance(n,swt) == 0))
              { # Get the true variance {var}:
                var = weight_variance(n,swt);
                if (abs(var - variance) < abs(bvar - variance))
                  { bvar = var; bxvr = xvr; bfac = fac; }
                n_accept++;
              }
            else
              { n_reject++; }
          }
      }
    printf "rejected %d  accepted %d\n", n_reject, n_accept > "/dev/stderr";
    # Recompute table with optimal parameters {var,fac}: 
    printf "optimal  xvr = %10.8f fac = %6.2f\n", bxvr, bfac > "/dev/stderr";
    init_weights(bxvr,n,owt);
    normalize_weights(n,bfac,owt,swt,balanced);
    print_weights(n,swt);
  }

function abs(x)
  { return (x < 0 ? -x : x); }

function init_weights(var,n,wt,   k,xlo,xhi,sg)
  { # Fills {wt[0..n-1]} with gaussian-like weights of variance {var}.
    for (k = 0; 2*k < n; k++)
      { xlo = k - (n/2.0);
        xhi = xlo + 1;
        wt[k] = gaussian_integral(var,xlo,xhi);
        wt[n-1-k] = wt[k];
        # printf "[%s _ %s] gauss = %s \n", xlo, xhi, wt[k] > "/dev/stderr";
      }
  }

function gaussian_integral(var,xlo,xhi,   N,k,r,x,d2,sum)
  { # Integral of a Gaussian with deviation {variance} from {xlo} to {xhi}.
    N = 5; # Number of samples to take
    sum = 0;
    for (k = 0; k < N; k++)
      { r = (k + 0.5)/N;
        x = (1 - r)*xlo + r*xhi;
        d2 = (x*x)/var/2;
        if (d2 < 50) { sum += exp(-d2); }
      }
    return sum/N;
  }

function normalize_weights(n,fac,owt,swt,balanced,   k,i,sum,v,big)
  { # Saves in {swt} a variant on the weights {owt}, balanced if needed,
    # then scaled so that the largest element is {fac},
    # and converted to integers. 
    
    # If {balanced} is true, the weights are adjusted so that
    # the sum of even elements is equal to the sum of odd elements. 
    # Then {swt} is a partition-of-unity when replicated with stride 2.
    
    # Copy the weights:
    for(k = 0; k < n; k++) { swt[k] = owt[k]; }
    
    if (balanced)
      { # Compute {sum[r]} = sum of elements with parity {r}
        split("",sum); sum[0] = 0; sum[1] = 0;
        for(k = 0; k < n; k++) { sum[k%2] += swt[k]; }

        # Normalize each element so that the 2-stride sums are 1:
        for(k = 0; k < n; k++) { swt[k] /= sum[k%2]; }
      }

    # Find largest element {big}:
    big = 0;
    for(k = 0; k < n; k++) 
      { v = swt[k];
        if (v < 0) { v = -v; }
        if (v > big) { big = v; }
      }
    
    # Convert to ints, scaled by {fac/big}:
    for(k = 0; k < n; k++) 
      { swt[k] = int(swt[k]*fac/big + 0.5);
        if ((k >= n-1-k) && (swt[k] != swt[n-1-k]))
          { printf "oops(3) swt[%d] = %s  !=  swt[%d] = %s -- owt[%d] = %s owt[%d] = %s\n", \
              k, swt[k], n-1-k, swt[n-1-k], \
              k, owt[k], n-1-k, owt[n-1-k] \
              > "/dev/stderr"; 
          }
      }
  }

function weight_unbalance(n,wt,  k,sum)
  { # Computes the difference between
    # the odd-indexed and even-indexed weights. 
    split("",sum); sum[0] = 0; sum[1] = 0;
    for(k = 0; k < n; k++) { sum[k%2] += wt[k]; }
    # printf "sum[0] = %s sum[1] = %s\n", sum[0], sum[1] > "/dev/stderr";
    return sum[1] - sum[0];
  }

function weight_variance(n,wt,  k,sum,avg,d,var)
  { # Computes the actual variance of the weights {wt[0..n-1]}.
    # Compute the mean {avg}:
    sum = 0.0; avg = 0.0;
    for(k = 0; k < n; k++) { sum += wt[k]; avg += wt[k]*k; }
    avg = avg/sum; 
    # Compute the variance {var}:
    var = 0.0;
    for(k = 0; k < n; k++) { d = k - avg; var += wt[k]*d*d; }
    var = var/sum;
    # Return the standard deviation:
    return var;
  }

function print_weights(n,wt,   k,rin,ref,skp)
  { # Prints the weights {wt[0..n-1]} on one line, with comment showing radius and variance.
    printf "# Last edited on DATE TIME by USER\n";
    # Cmpute the input radius {rin} and effective radius {ref}:
    rin = (n-1)/2;
    skp = 0;
    while ((skp < rin) && (wt[skp] == 0) && (wt[n-1-skp] == 0)) { skp++; }
    ref = rin - skp;
    # Print the effective radius and variance:
    printf "# r = %2d  var = %7.4f\n", ref, weight_variance(n,wt);
    printf "\n";
    # Printe the weights, excluding zero tails:
    for(k = skp; k < n-skp; k++) { printf " %4d", wt[k]; }
    printf "\n";
  }

function arg_error(msg)
  { printf "** %s\n", msg;
    abort = 1;
    exit(abort);
  }
