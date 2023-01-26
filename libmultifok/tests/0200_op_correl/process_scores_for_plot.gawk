#! /usr/bin/gawk -f
# Last edited on 2023-01-25 08:00:26 by stolfi

# Reads "-cdata.txt" file produced by {mf_0200_op_correl.c}
# Writes a file with lines "{ki} {zave} {sharp} {score/scoreMax} {error}" for plotting.
# where {ki} is the image index, 
# the {sharp} and {score} are averages binned by {Z},
# and {error} is {score-sharp}.
# The {scoreMax} is the max score in each image.
# 
# Blank lines are written between images.

BEGIN{ 
  abort = -1;
  split("", sum_wp_sharp); # Sum of {sharp*wp} per histogram bin.
  split("", sum_wp_score); # Sum of {score*wp} per histogram bin.
  split("", sum_wp);       # Sum of {wp} per histogram bin.
  split("", sharps);       # Histogram values of {sharp}.
  split("", scores);       # Histogram values of {score}.
  
  # Parameters of histogram:
  zaveMin = -30;
  zaveMax = +30;
  zaveStep = 1.00;
  nh = ceil((zaveMax - zaveMin)/zaveStep) + 1; # Max number of hist bins.
  # Curreny image number:
  oki = -1;
  
}

(abort >= 0) { exit(abort); }

/^ *[P][0-9]/ {
  id = $1;    # Image and pixel id.
  wp = $2;    # Pixel weight for regression.
  sharp = $3; # Actual sharpness.
  score = $4; # Computed sharpness score.
  zave = $5;  # Average {Z} coord in pixel relative to focus plane height.
  zdev = $6;  # Deviation of {Z} coord in pixel.
  
  # Extract the image index {ki}: 
  ki = id; gsub(/[P]/,"",ki); gsub(/[.].*/,"",ki); ki += 0;
  if (ki != oki) {
    # Starting data for a new image.
    if (oki != -1) { write_image_data(oki); }
    clear_image_data(ki);
    oki = ki;
  }
  accumulate_image_data(ki, wp, sharp, score, zave, zdev);
  next;
}

// { data_error("invalid line format"); }

END {
  if (abort >= 0) { exit(abort); }
  if (oki != -1) { write_image_data(oki); }
}

function clear_image_data(ki,   ih) {
  for (ih = 0; ih < nh; ih++) {
    sum_wp_sharp[ih] = 0;
    sum_wp_score[ih] = 0;
    sum_wp[ih] = 0;
    np = 0; # Number of data lines read.
    nu = 0; # Number of data lines used in plot.
  }
}
  
function accumulate_image_data(ki,wp,sharp,score,zave,zdev, ih) {
  np++;
  if (zdev < 0.3) {
    ih = floor((zave - zaveMin)/zaveStep);
    sum_wp_sharp[ih] += wp*sharp;
    sum_wp_score[ih] += wp*score;
    sum_wp[ih] += wp;
    nu++;
  }
}

function write_image_data(ki, ih,scoreMax,scoreMin,zaveih,nw) {
  scoreMax = -9999;
  scoreMin = 9999;
  for (ih = 0; ih < nh; ih++) {
    if (sum_wp[ih] > 0) {
      sharps[ih] = sum_wp_sharp[ih]/sum_wp[ih];
      scores[ih] = sum_wp_score[ih]/sum_wp[ih];
      if(scores[ih] > scoreMax) { scoreMax = scores[ih]; }
      if(scores[ih] < scoreMin) { scoreMin = scores[ih]; }
    }
  }
  nw = 0;
  for (ih = 0; ih < nh; ih++) {
    if (sum_wp[ih] > 0) {
      zaveih = zaveMin + (ih + 0.5)*zaveStep;
      printf "%+12.6f", zaveih; 
      printf "  %12.6f %+12.6f", sharps[ih], scores[ih]/scoreMax; 
      printf "  %+12.6f\n", scores[ih] - sharps[ih];
      nw++;
    } else {
      printf "\n"; # Break in plot.
    }
  }
  printf "\n";
  printf "image %d: %d lines read, %d used, %d discarded", ki, np, nu, np-nu > "/dev/stderr"; 
  printf " %d plot points", nw > "/dev/stderr"; 
  printf " -- score range [ %+12.6f _ %+12.6f ]", scoreMin, scoreMax > "/dev/stderr"; 
  printf " ratio %+12.6f\n", scoreMin/scoreMax > "/dev/stderr"; 
}
 
function floor(x)
  { x = (x >= 0? int(x) : -int(-x));
    return x;
  }
 
function ceil(x)
  { x = (x <= 0? int(x) : -int(-x));
    return x;
  }
         
function data_error(msg)
  { printf "%s:%s: ** %s\n", FILENAME, FNR, msg > "/dev/stderr"; 
    printf "  «%s»\n", $0 > "/dev/stderr"; 
    abort = 1;
    exit(abort);
  } 
          
function arg_error(msg)
  { printf "** %s\n", msg > "/dev/stderr"; 
    abort = 1;
    exit(abort);
  } 

