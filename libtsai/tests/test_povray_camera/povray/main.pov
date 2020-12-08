#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"
#include "luzes.inc"



#macro ponto(R, tx, ty, tz)
  union{
    sphere{
      <0,0,0>,R
      texture{
        pigment{ color rgb <1.000,1.000,1.000> }
        finish{ ambient 0.50 diffuse 0.50 }
      }
    }
    translate <tx,ty,tz>
  }
#end



#macro axis(R,r,dir,cor)
  union{
    cylinder{
      <0,0,0>,R*dir,r
     
    }
    #local step = 50;
    #local i = step;
    #while (i < R) 
        sphere{ i*dir,(3+2*mod(i/step,2))*r}
        #local i = i+step;
    #end 
    texture{
        pigment{ color rgb cor }
        finish{ ambient 0.50 diffuse 0.50 }
    }    
  }
#end


#macro axes(R,r)
  union{
    axis(R,r,x,<1,0,0>)
    axis(R,r,y,<0,1,0>)
    axis(R,r,z,<0,0,1>) 
    ponto(6*r, 0,0,0)
  }
#end

// Determina parâmetros da câmera em funcao do tempo
#declare ctr = <0,0,0>;  // Centro de interesse

luzes(1.0)
  
axes(3000,3)


// Define a posicao da camera:
#include concat("camera_dat.inc")
