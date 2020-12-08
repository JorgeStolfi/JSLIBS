#! /bin/bash
# Last edited on 2012-01-10 02:04:54 by stolfilocal

na="$1"; shift;
dfile="$1"; shift;
d0="0.25"
d1="0.50"
d2="1.00"
d3="2.00"
d4="4.00"
d5="8.00"

gnuplot << EOF
  set term x11
  awd(d) = (1/d)/sqrt((0.50/d)**2 + 1)
  # set xrange [0.1:10.0]
  # set yrange [-0.001: 10.00]
  # plot (awd(x)) title "awd(d)" with lines
  # pause 30
  na = ${na}
  awh(d) = sqrt((awd(d))**2 + (0.5*pi/na)**2)
  gau(x,y,sigma) = exp(-0.5*((x-1)**2 + y**2)/sigma**2)/(sigma*sqrt(pi))
  ref(a,d) = (4.35/na)*gau(cos(a),sin(a),awh(d))
 
  # distance from {(d,0)} to {(r*cos(a),r*sin(a))}:
  dst(a,d,r,m) = sqrt((m+(r-m)*cos(a)-d)**2 + (r*sin(a))**2)
  
  # distance squared from {(1,0)} to {(cos(a),sin(a))} with corrections:
  crr(a,m,n) = m/4*(1-cos(a))**2 - n/8*(1-cos(a))**3
  ds2(a,d,m,n) = ((sin(a))**2) + (1-cos(a))**2 - crr(a,m,n)

  # model effective distance squared for distance {d} with correction coeffs {m,n}:
  crr(a,m,n) = m/4*(1-cos(a))**2 - n/8*(1-cos(a))**3
  def(a,d,g,h) = ds2(a,d,0,0) - g*(ds2(a,d,0,0))**2 + h*(ds2(a,d,0,0))**3

  # effective distance squared computed from probability:
  dpr(id,e,p) = -log(column(3+id)/e)*2*(p**2)
  
  d0 = ${d0}
  d1 = ${d1}
  d2 = ${d2}
  d3 = ${d3}
  d4 = ${d4}
  d5 = ${d5}
  
  # alpha = 0.0174532922


  e0 =  1.3405739012
  e1 =  1.7436732929
  e2 =  2.7044742723
  e3 =  4.9994280400
  e4 =  9.8888214319
  e5 = 19.5096541885
  
  p0 = 1.7111
  p1 = 1.1583
  p2 = 0.7502
  p3 = 0.4480
  p4 = 0.2425
  p5 = 0.1240
  
  m0 = 0 # 0.347
  m1 = 0 # 0.697
  m2 = 0 # 1.390
  m3 = 0 # 2.574
  m4 = 0 # 3.770
  m5 = 0 # 4.432

  n0 = 0 # 0.022
  n1 = 0 # 0.080
  n2 = 0 # 0.279
  n3 = 0 # 0.772
  n4 = 0 # 1.330
  n5 = 0 # 1.632
  
  g0 = 0
  g1 = 0
  g2 = 0
  g3 = 0
  g4 = 0
  g5 = 0

  h0 = 0 
  h1 = 0 
  h2 = 0 
  h3 = 0 
  h4 = 0 
  h5 = 0

  dfile="${dfile}"

  call "plot-polar-gaussian-1.gpl"
  pause 120

  call "plot-polar-gaussian-3.gpl"
  pause 120

  call "plot-polar-gaussian-2.gpl"
  pause 120

  quit
EOF
