// Last edited on 2010-05-01 19:03:48 by stolfilocal

global_settings { assumed_gamma 2.2 }

#declare closeup = false; // Shows close-up view if true; full view if false.

background { color rgb <1.000, 1.000, 1.000>  }
#include "scene_Z.inc"

// Radius (half-diagonal) of domain rectangle:
#declare radiusXY = sqrt(maxX*maxX + maxY*maxY)/2;

#declare tx_plot_color = < 0.300, 1.000, 0.500 >;

#declare tx_plot = 
  texture{ 
    pigment{ color rgb tx_plot_color }  
    finish {diffuse 0.8 ambient 0.2} 
  }

union{
  object{ scene_mesh translate < -(maxX)/2.0,-(maxY)/2.0,0 > scale < 1, 1, scaleZ > } 
  // sphere{ <0,0,0>, 0.5*radiusXY }
  texture { tx_plot }
}

#if (closeup)
  // For debugging. Edit these as needed.
  #declare viewAzim = 110.0;
  #declare viewElev = 30.0;
  #declare camera_dir = vrotate(vrotate(<1,0,0>, -viewElev*y), -viewAzim*z);
  #declare camera_disp = 1000*camera_dir;
  #declare camera_radius = 3.2;
  #declare camera_focus = <26.0 - (maxX)/2.0, 118.0 - (maxY)/2.0, -30.552839*scaleZ >;
#else
  // For production figures. Edit with care.
  #declare camera_dir = vrotate(vrotate(<1,0,0>, -viewElev*y), -viewAzim*z);
  #declare camera_disp = viewDist*camera_dir;
  #declare camera_radius = 1.15*radiusXY;
  #declare camera_focus = <0,0,focusZ*scaleZ>;
#end
#declare image_rel_width = image_width;
#declare image_rel_height = image_height;

#include "camera.inc"

union{
  light_source { 
    <20000, 0, 0> color rgb 0.60*<1,1,1> 
    rotate -(+50 + viewElev/2)*y
    rotate -(+30 + viewAzim)*z
  }
  light_source { 
    <20000, 0, 0> color rgb 0.25*<1,1,1>
    rotate -(+10 + viewElev/2)*y
    rotate -(+05 + viewAzim)*z
  }
  light_source { 
    <20000, 0, 0> color rgb 0.25*<1,1,1>
    rotate -(+22 + viewElev/2)*y
    rotate -(-20 + viewAzim)*z
  }
}
