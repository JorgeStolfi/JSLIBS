#! /usr/bin/gawk -f
# Last edited on 2023-01-30 12:18:50 by stolfi

# Determines the terms for quadratic regression with the {CANC} or {HART} basis coefficient
# for a 3x3 window.
#
# The user must define (with "-v") the variable {belNames} which is the names
# of the {nb=9} elements of the basis, in the same order as used by {multifok_basis_make},
# separated by spaces.  The elements must be named "{L}{x}{y}" where {L} is "X" or "H",
# and {x} and {y} are "m", "o", or "p" meaning {-1}, 0, or {+1}, respectively.
#
# The user must also define the {outprefix} for output file names.
#
# The program writes a term names file and a term indices file.
# 
# The term names file is "{outPrefix}-terms.txt". Each line has a term
# index {kt} from 0 to {nt-1}, and a sum of one or more products
# "{el1}*{el2}" of pairs of basis element names, such as "Xoo*Xmo +
# Xoo*Xom + Xoo*Xop + Xoo*Xpo".
#
# The products that appear are only those with are distinct apart from
# communativity; meaning that only one of "Xop*Xmm" and "Xmm*Xop" will
# appear. Thus there are only {np-nb*(nb+1)/2=45} distinct products.
# Each product is listed in its canonical form, namely with {el1} not
# after {el2} in the given order of basis elements.
#
# Each term has all and only the canonical products that are equivalent under
# axis flips ("m" <--> "p") and transposition of indices ("mo" <--> "om", etc.)
# followed by re-canonicalization of the product.  Each canonical product
# appears in exaclty one term. The number of terms {nt} is determined by the program.
#
# The term indices file is "{outPrefix}-tix.txt". Each line corresponds
# to a product "{el1}*{el2}". It has the format "{ip} {jb1} {jb2} {kt}
# {pr}" where {ip} is the seqeuntial product index from 0 to {np-1},
# {jb1} and {jb1} are the indices of the two basis elements whose
# coefficients are to be multiplied, {kt} is the index of the term to
# which the product is to be added, and {pr} is the textual product,
# e.g. "Xmo*Xmp".


BEGIN { 
  if (belNames == "") { arg_error("must define {belNames}"); }
  if (outPrefix == "") { arg_error("must define {outPrefix}"); }

  # Define the list {el_jb[1..nb]} of basis element names in canonical (radial) order.
  # Note indices start at 1.
  nb = split(belNames, el_jb); 
  if (nb != 9) { arg_error(("invalid basis size = " nb)); }
  split("", jb_el); # Maps element name to index in {el_jb}.
  for (jb = 1; jb <= nb; jb++) { 
    el = el_jb[jb];
    if (el !~ /^[SH][mop][mop]$/)  { arg_error(("invalid basis elem name = " el)); }
    jb_el[el] = jb;
  }
  
  printf "generating prods...\n" > "/dev/stderr";
  np = nb*(nb+1)/2;    # Number of canonical products of pair of basis elem coeffs.
  split("", pr_ip);    # All canonical prods "{el1}*{el2}", indexed {0..np-1}.
  split("", jb1_ip);   # {el_jb[jb1[ip]]} is the first element of the product.
  split("", jb2_ip);   # {el_jb[jb2[ip]]} is the second element of the product.
  split("", ip_pr);    # Maps canonical prod name {pr} to index in {pr_ip}.
  get_all_canonical_products(); 
  
  printf "assigning products to terms...\n" > "/dev/stderr";
  # For {kt} in {0..nt-1}, {nr_kt[kt]} is count of products of term {kt}. The indices of those products
  # are {ip_kt_r[kt,r]} for {r} in {0..nr_kt[kt]-1}.
  #
  nt = 0; # Number of terms created. 
  split("", nr_kt);    # 
  split("", ip_kt_r); # 
  split("", kt_ip);    # Maps canonical prod index {ip} to the index {kt} of the term that includes it.
  nt = assign_products_to_terms(); 
    
  printf "obtained %d terms\n", nt > "/dev/stderr";
  write_term_names_file();
  write_term_indices_file();
}

function write_term_names_file(  fna,kt,nr,r,ip,pr,sep) {
  fna = (outPrefix "-terms.txt");
  for (kt = 0; kt < nt; kt++) {
    printf "%3d ", kt > fna;
    sep=""
    if (! (kt in nr_kt)) { prog_error("bug {nr_kt[kt]} undef"); }
    nr = nr_kt[kt]; # Number of products of this term.
    if (nr < 0) { prog_error("bug {nr_kt[kt]} invalid"); }
    for (r = 0; r < nr; r++) { 
      ip = ip_kt_r[kt,r];
      pr = pr_ip[ip];
      printf "%s%s", sep, pr > fna;
      sep="+";
    }
    printf "\n" > fna;
  }
  close(fna);
}
   
function write_term_indices_file(  fna,ip,jb1,jb2,pr,kt) {
  fna=(outPrefix "-tix.txt");
  for (ip = 0; ip < np; ip++) {
    jb1 = jb1_ip[ip];
    jb2 = jb2_ip[ip];
    kt = kt_ip[ip]
    pr = pr_ip[ip];
    printf "%3d  %3d %3d  %3d %s\n", ip, jb1-1, jb2-1, kt, pr > fna;
  }
  close(fna);
}

function get_all_canonical_products(   jb1,dj,jb2,el1,el2,ip,pr,red) {
  # Generates all canonical products, starting with the squares.
  ip = 0;
  for (dj = 0; dj < nb; dj++) {
    for (jb2 = 1 + dj; jb2 <= nb; jb2++) {
      jb1 = jb2 - dj;
      if ((jb1 < 1) || (jb1 > jb2) || (jb2 > nb)) { prog_error("bug {jb1,jb2} enum"); }
      # Another canonical product:
      el1 = el_jb[jb1];
      el2 = el_jb[jb2];
      pr = (el1 "*" el2)
      # Store in the product table: 
      pr_ip[ip] = pr;
      jb1_ip[ip] = jb1;
      jb2_ip[ip] = jb2;
      ip_pr[pr] = ip;
      printf "  product %3d = (%d,%d) = %s\n", ip, jb1,jb2, pr > "/dev/stderr";
      ip++;
    }
  }
  if (ip != np) { prog_error("bug {kp}"); }
}

function assign_products_to_terms(  ip,pr,jb1,jb2,pr_red,ip_red,kt_red,nt) {
  nt = 0;
  # For safety: 
  for (ip = 0; ip < np; ip++) { kt_ip[ip] = -1; }
  # Now find equivalences:
  for (ip = 0; ip < np; ip++) { 
    pr = pr_ip[ip];
    jb1 = jb1_ip[ip];
    jb2 = jb2_ip[ip];
    printf "product %3d = %d,%d = %s", ip, jb1,jb2, pr > "/dev/stderr";
    ip_red = least_equivalent_prod(jb1,jb2);
    pr_red = pr_ip[ip_red];
    printf " --> %3d = %s\n", ip_red, pr_red > "/dev/stderr";
    if (kt_ip[ip_red] < 0) {
      # This is the lowest equiv product of a new term.
      # Assign product {ip_red} to a new term:
      kt = nt; # Index of the new term.
      r = 0;   # Index of the product {ip_red} in that term.
      kt_ip[ip_red] = kt;
      ip_kt_r[kt,r] = ip_red;
      nr_kt[kt] = 1;
      nt++;
    } 
    if (ip != ip_red) {
      # Assign product {ip} to the same term as {ip_red}:
      kt_red = kt_ip[ip_red];
      kt_ip[ip] = kt_red;
      r = nr_kt[kt_red];
      ip_kt_r[kt_red,r] = ip;
      nr_kt[kt_red]++;
    } else {
      if (kt_ip[ip] < 0) { prog_error("bug {kt_ip[ip]}"); }
    }
  }
  return nt;
}

function least_equivalent_prod(jb1,jb2,   ip_red,ip,ix) {
  # for (ix in ip_pr) { printf "    ip_pr[%s] = %s\n", ix, ip_pr[ix] > "/dev/stderr"; }
  ip_red = 9999;
  # We enumerate all pairs {jb1,jb2} equivalent to the given pair,
  # and choose the one which gives a product with lowest {ip}.
  for (ks = 0; ks < 2; ks++) {
    # Transpose loop.
    for (kv = 0; kv < 2; kv++) {
      # Vertical flip loop.
      for (kh = 0; kh < 2; kh++) {
        # Horizontal flip loop.
        # At this point "{el_jb[jb1]}*{el_jb[jb2]}" is a product eqiv to the given one.
        # Get its index {ip} accounting for commutativity canonicalization:
        if (jb1 <= jb2) {
          ip = ip_from_jb1_jb2(jb1,jb2)
        } else {
          ip = ip_from_jb1_jb2(jb2,jb1)
        }
        if ((ip < ip_red)) { ip_red = ip; }
        # Apply a horizontal flip to the image:
        jb1 = flip_hor(jb1);
        jb2 = flip_hor(jb2);
      }
      # Apply a vertical flip to the image:
      jb1 = flip_ver(jb1);
      jb2 = flip_ver(jb2);
    }
    # Apply an axis swap to the image
    jb1 = axis_swap(jb1);
    jb2 = axis_swap(jb2);
  }
  return ip_red;
}

function ip_from_jb1_jb2(jb1,jb2,  el1,el2,pr,ip,ix) {
  # for (ix in ip_pr) { printf "    ip_pr[%s] = %s\n", ix, ip_pr[ix] > "/dev/stderr"; }
  el1 = el_jb[jb1];
  el2 = el_jb[jb2];
  pr = (el1 "*" el2);
  if (! (pr in ip_pr)) { prog_error(("invalid product " jb1 "," jb2 " = " pr)); }
  return ip_pr[pr];
}

function flip_hor(jb, el,L,x,y) {
  el = el_jb[jb];
  L = substr(el,1,1);
  x = substr(el,2,1);
  y = substr(el,3,1);
  el = (L (x == "m" ? "p" : (x == "p" ? "m" : x)) y);
  return jb_el[el];
}

function flip_ver(jb, L,x,y) {
  el = el_jb[jb];
  L = substr(el,1,1);
  x = substr(el,2,1);
  y = substr(el,3,1);
  el = (L x (y == "m" ? "p" : (y == "p" ? "m" : y)));
  return jb_el[el];
}

function axis_swap(jb, L,x,y) {
  el = el_jb[jb];
  L = substr(el,1,1);
  x = substr(el,2,1);
  y = substr(el,3,1);
  el = (L y x);
  return jb_el[el];
}
         
function data_error(msg) { 
  printf "%s:%s: ** %s\n", FILENAME, FNR, msg > "/dev/stderr"; 
  printf "  «%s»\n", $0 > "/dev/stderr"; 
  abort = 1;
  exit(abort);
} 
          
function arg_error(msg) { 
  printf "** %s\n", msg > "/dev/stderr"; 
  abort = 1;
  exit(abort);
} 
          
function prog_error(msg) { 
  printf "** prog error: %s\n", msg > "/dev/stderr"; 
  abort = 1;
  exit(abort);
} 
