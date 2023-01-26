#! /usr/bin/gawk -f
# Last edited on 2023-01-24 19:46:15 by stolfi

BEGIN {
  split("", sc_ds);
  split("", sc_wp);
  np = 0;
}

/^ *[P]/ {
  id = $1; wp = $2; sharp = $3; score = $4;
  sc_ds[np] = score-sharp;
  sc_wp[np] = wp;
  np++;
  next;
}

// { 
  printf "** line %d: bad format [[%s]]\n", FNR, $0 > "/dev/stderr"; 
  exit(1);
}

END {
  sum_wp = 1.0e-200;
  sum_wp_ds = 0;
  for (kp = 0; kp < np; kp++)
    { ds = sc_ds[kp];
      wp = sc_wp[kp];
      sum_wp_ds += wp*ds;
      sum_wp += wp;
    }
  avg = sum_wp_ds/sum_wp;
  sum_wp_ds2 = 0;
  for (kp = 0; kp < np; kp++)
    { ds = sc_ds[kp];
      wp = sc_wp[kp];
      sum_wp_ds2 += wp*ds*ds;
    }
  dev = sqrt(sum_wp_ds2/sum_wp);
  printf "AVG error {score - sharp} = %12.4f\n", avg > "/dev/stderr";
  printf "RMS error {score - sharp - AVG} = %12.4f\n", dev > "/dev/stderr";
}
