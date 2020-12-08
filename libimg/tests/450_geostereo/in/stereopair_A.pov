// Last edited on 2017-06-26 04:54:42 by stolfilocal
// Processed by remove-cam-lights

background{ color rgb < 1.00, 0.60, 0.85 > }

#include "params.inc"

#declare tx_plastico =
  texture{
    pigment{ color rgb < 0.11, 0.80, 1.00 > }
    finish{ diffuse 0.8 ambient 0.1 specular 0.5 roughness 0.005 }
  }

#declare tx_fosca =
  texture{
    pigment{ color rgb < 1.00, 0.80, 0.10 > }
    finish{ diffuse 0.9 ambient 0.1 }
  }

#declare tx_cinza =
  texture{
    pigment{ color rgb < 0.75, 0.75, 0.75 > }
    finish{ diffuse 0.9 ambient 0.1 }
  }

#declare tx_xadrez =
  texture{
    pigment{ checker color rgb < 0.10, 0.32, 0.60 >, color rgb < 1.00, 0.97, 0.90 > }
    finish{ diffuse 0.9 ambient 0.1 }
    scale 2.0
  }

#macro interpola_1(tt,tt0,vv0,tt1,vv1)
  #local c = ( (tt-tt0)/(tt1-tt0) );
  (vv0 + c*(vv1-vv0))
#end

#declare die = seed(4615);

#macro randpt(rad)
  #local len = 10;
  #while (len > 1.0)
    #local pt = 2 * < rand(die), rand(die), rand(die) > - <1,1,1>;
    #local len = vlength(pt);
  #end
  rad*pt
#end

#macro gera_palitos(N,rtangle,rr)
  #local amp = 0.3;
  // #local pi = 3.1415926;;
  union{
    #local k = 0;
    #while (k <= N)
      #local kf = k/N; // Phase [0 - 1]
      #local x0 = 2*(k + rand(die))/N - 1.0;
      #local x1 = 2*(k + rand(die))/N - 1.0;
      #local p0 = rtangle * < x0, -1.0, amp*sin(2*pi*kf) >;
      #local p1 = rtangle * < x1, +1.0, amp*sin(3*2*pi*kf) >;
      cylinder{ p0, p1, rr texture{ tx_plastico } }
      #local k = k + 1;
    #end
  }
#end

#declare chao =
  box{ <-20,-20,-1>, <+20,+20,0> }

#include "eixos.inc"

#declare rtangle = 150;

union{
  object{ gera_palitos(200,150,0.75) }
  // object { eixos(150.0) }
  // object {chao    texture{tx_xadrez}}
}

#include "camlight.inc"
#declare dir_camera = < 2.00, -3.00, 14.00 >;
#declare dir_dir = vnormalize(< -dir_camera.y, dir_camera.x, 0.0>);
#declare centro_cena = < 0.00, 0.00, 0.00 > + 3.5*(2*eye - 1)*dir_dir;
#declare raio_cena = 2*rtangle;
#declare dist_camera = 10*raio_cena;
#declare intens_luz = 1.00;
camlight(centro_cena, raio_cena, dir_camera, dist_camera , z, intens_luz)
