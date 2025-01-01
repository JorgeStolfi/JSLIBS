// Last edited on 2024-11-08 17:22:54 by stolfi

#include "lamp_array.inc"
#include "camlight.inc"

background{ color rgb < 0.900, 0.970, 0.990 > }

#macro make_texture(C)
  texture{ 
    pigment{ color rgb 0.0*C }
    finish{ diffuse 0.60 ambient 0.30 }
  }
  texture{ 
    pigment{ color rgb <0,0,0> filter 1.0 }
    finish{ diffuse 0 ambient 0 reflection 0.30*C }
    // normal{ granite 0.5 scale 0.1 }
  }
#end

#declare tx_amarelo = make_texture(< 0.980, 1.000, 0.200 >)
#declare tx_vermelho = make_texture(< 1.000, 0.250, 0.200 >)
#declare tx_branco = make_texture(< 1.000, 1.000, 1.000 >)
  
#macro electron(P)
  sphere{ P, 0.08 texture{ tx_amarelo } } 
#end
    
#macro proton(P)
  sphere{ P, 0.10 texture{ tx_vermelho } } 
#end

#macro platter()
  cylinder{ -0.1*z, 0*z, 2.0 translate -1.25*z texture{ tx_branco } }
#end

union{
  object{ platter() }
  object{ electron( <0.8, 0.8, 0.2> ) }  
  object{ proton( <0,0,0> ) }
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
