#! /usr/bin/gawk -f
# Last edited on 2013-10-11 16:38:07 by stolfilocal

BEGIN \
  { abort = -1;
    usage = ( "compute-weights-variance \\\n" \
      "  < INFILE \\\n" \
      "  > OUTFILE" \
    );

    # For each line that starts with a "{XXX}WEIGHTS := " (possibly
    # preceded by "#") collects the numbers after the ":=" and
    # computes their variance {VAR} and even/odd balance {UNB}. Then
    # appends the comment "# var ={VAR} unb = {UNB}" at the end of the
    # line. If the comment is already there, simply updates the values
    # of {VAR} and {UNB}.
}

/^[\#]?[ ]*[A-Za-z_]+WEIGHTS[ ][:]?[=]/ {
  # Weights definition line.
  # Save the line:
  lin = $0;
  # Remove everything but the weights:
  gsub(/^[\#]?[ ]*[A-Za-z_]+WEIGHTS[ ][:]?[=]/, "", $0);
  gsub(/[ ]*[\#].*$/, "", $0);
  # Grab the weights:
  split("", wt);
  n = NF;
  for (i = 0; i < n; i++) { wt[i] = $(i+1); }
  # Compute their variance and unbalance:
  var = weight_variance(n,wt);
  unb = weight_unbalance(n,wt);
  # Add values to line:
  lin = append_value(lin, "var =", sprintf("%6.4f", var));
  lin = append_value(lin, "unb =", sprintf("%+6.4f", unb));
  # Output line:
  print lin;
  next;
}

// { 
  # Other line; just keep unchanged.
  print; next;
}

function abs(x)
  { return (x < 0 ? -x : x); }

function append_value(lin,tag,val,   cmt,i,ini,fin)
  { i = index(lin,tag);
    if (i > 0)
      { # The {tag} is present in {lin}; update its value
        ini = substr(lin,1,i-1);         # Part of {lin} up to the {tag}
        fin = substr(lin,i+length(tag)); # Part of {lin} after the {tag}
        gsub(/^[ ]*[-+]?[0-9]*([0-9]|[.])[0-9]*/, "", fin);
        lin = (ini tag " " val fin);
      }
    else
      { # The {tag} is not present in {lin}.
        if (! match(lin, /[^ ].*[\#]/))
          { # Line {lin} does not have a '#' except possibly at beginning.
            lin = (lin "#");
          }
        lin = (lin " " tag " " val);
      }
    return lin;
  }

function weight_unbalance(n,wt,  k,sum)
  { # Computes the difference between
    # the odd-indexed and even-indexed weights. 
    split("",sum); sum[0] = 0; sum[1] = 0;
    for(k = 0; k < n; k++) { sum[k%2] += wt[k]; }
    # printf "sum[0] = %s sum[1] = %s\n", sum[0], sum[1] > "/dev/stderr";
    return (sum[1] - sum[0])/(sum[1] + sum[0]);
  }

function weight_variance(n,wt, k,sum,avg,d,var)
  { # Computes the actual variance of the weights {wt[0..n-1]}
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

function print_weights(n,wt,   k)
  { # Print the weights, scaled by 999 and rounded to ints:
    printf "# F_WEIGHTS :=  ";
    for(k = 0; k < n; k++) { printf " %3d", wt[k]; }
    # Print the standard deviation:
    printf "  # sigma = %6.4f K = \n", weight_deviation(n,wt);
  }

function arg_error(msg)
  { printf "** %s\n", msg;
    abort = 1;
    exit(abort);
  }
