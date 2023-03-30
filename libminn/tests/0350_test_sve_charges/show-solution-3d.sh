#! /bin/bash
# Last edited on 2012-01-22 23:00:16 by stolfilocal

# Plots the electron positions found by {test_sve_charges}

nq="$1"; shift;
Lmax="$1"; shift;
solfile="$1"; shift;

# Remove leading zeros from ${nq} (THE OCTAL CONVENTION IS A CROCK!!!!!!!)
nq=`echo "${nq}" | sed -e 's:^0*::'`

povfile=${solfile%%.*}.pov
imgfile=${solfile%%.*}.png

if [[ -e /usr/bin/povray ]]; then
  POVRAY="/usr/bin/povray"
elif [[ -e /usr/local/bin/povray ]]; then
  POVRAY="/usr/local/bin/povray"
else
  echo "no povray?" 1>&2; exit 1
fi
POVINC1="/usr/share/povray/include"
POVINC2="povray"
POVTTF="povray/ttf"

cat > ${povfile} <<EOF

#include "lamp_array.inc"
#include "camlight.inc"

background{ color rgb < 0.900, 0.970, 0.990 > }

#declare tx_cloud = 
  texture{ 
    pigment{ color rgb < 0.950, 0.980, 0.990 > filter 0.900 }
    finish{ diffuse 0.00 ambient 0.00 specular 0.20 roughness 0.05 }
  }

#declare tx_electron = 
  texture{ 
    pigment{ color rgb < 0.980, 1.000, 0.200 > }
    finish{ diffuse 0.60 ambient 0.30 specular 0.10 roughness 0.005 }
  }

#declare tx_bar = 
  texture{ 
    pigment{ color rgb < 0.200, 0.700, 0.200 > }
    finish{ diffuse 0.60 ambient 0.30 specular 0.10 roughness 0.005 }
  }

#declare tx_platter = 
  texture{ 
    pigment{ color rgb < 1.000, 1.000, 1.000 > }
    finish{ diffuse 0.90 ambient 0.10 }
  }
  
#declare cloud_R =    1.00;
#declare electron_R = 0.06;
#declare bar_R =      0.01;
#declare bar_Lmax =   ${Lmax};
#declare platter_R =  1.25;
#declare platter_th = 1.25;
   
#macro cloud()
  sphere{ <0,0,0>, cloud_R texture{ tx_cloud } } 
#end
 
#macro electron(P)
  sphere{ P, electron_R texture{ tx_electron } } 
#end
  
#macro bar(P,Q)
  #local L = vlength(P-Q);
  #if (L < bar_Lmax*sqrt(vlength(P)*vlength(P) + vlength(Q)*vlength(Q)))
    #local CH = mod(3*L,1);
    #local CR = 0.5*(1 + cos(2*pi*CH));
    #local CB = 0.5*(1 + sin(2*pi*CH));
    #local CG = (0.635 - 0.1*CB - 0.3*CR)/0.6;
    cylinder{ P, Q, bar_R texture{ tx_bar pigment{ color rgb < CR, CG, CB > } } } 
  #end
#end

#macro platter()
  cylinder{ -platter_th*z, 0*z, platter_R translate -1.25*cloud_R*z texture{ tx_platter } }
#end

EOF

cat ${solfile} \
  | while [[ 1 ]]; do
      read i cx cy cz dorg rest;
      if [[ "/${i}" == "/#" ]]; then continue; fi;
      if [[ "/${cx}" == "/" ]]; then break; fi;
      printf "#declare P%03d = < %+8.5f, %+8.5f, %+8.5f >;\n" "${i}" "${cx}" "${cy}" "${cz}">> ${povfile} ;
    done
      
cat >> ${povfile} <<EOF
union{
  object{ platter() }
  object{ cloud() }
EOF

for i in `seq 0 $(( ${nq} - 1 ))`; do
  printf "  object{ electron(P%03d) }" "${i}" >> ${povfile} ;
  for j in `seq 0 $(( ${i} - 1 ))`; do
    printf "  bar(P%03d,P%03d)" "${i}" "${j}" >> ${povfile} ;
  done
done

cat >> ${povfile} <<EOF
}

#declare cam_ctr = <0,0,0>;
#declare cam_dir = <5,4,3>;
#declare cam_dst = 5.0;
#declare cam_rad = 1.5;
#declare cam_upp = z;

#declare cam_lux = 1.0;
#declare cam_lnr = 1;
#declare cam_lar = 2;
#declare cam_lct = 0;

camlight(cam_ctr,cam_rad,cam_dir,cam_dst,cam_upp, cam_lux,cam_lnr,cam_lar,cam_lct)
EOF

rm -fv ${imgfile}
${POVRAY} \
    +FN +Q9 +MB1 \
    +W512 +H512 \
    +AM1 +A0.0 +R2 \
    +D +SP32 +EP4 \
    +L${POVINC1} \
    +L${POVINC2} \
    +L${POVTTF} \
    +I${povfile} \
    +O${imgfile} \
  2>&1 | ./povray/povray-output-filter.gawk
display ${imgfile}
