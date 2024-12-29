// Last edited on 2010-05-01 22:41:13 by stolfilocal

global_settings { assumed_gamma 2.2 }

#declare closeup = false; // Shows close-up view if true; full view if false.
#declare bicolor = true; // True to paint pos & neg with different hues.

background { color rgb <1.000, 1.000, 1.000>  }
#include "scene_eZ.inc"

// Radius (half-diagonal) of domain rectangle:
#declare radiusXY = sqrt(maxX*maxX + maxY*maxY)/2;

#declare tx_plot_color1 = < 1.000, 0.500, 0.300 >;
#declare tx_plot_color2 = < 1.000, 0.900, 0.600 >;

#if (bicolor)
  #declare tx_plot = 
    texture{ 
      pigment{ 
        image_map{ ppm "isolines.ppm" }
        rotate -90*x
        translate -0.5*x
        scale < 10*radiusXY, 1, 200*stepZ >
      }  
      finish { diffuse 0.8 ambient 0.2 } 
    }
#else
  #declare tx_plot = 
    texture{ 
      pigment{ 
        checker color rgb tx_plot_color1, color rgb tx_plot_color2
        translate < 0.5, 0.5, 0.0 >
        scale < 10*radiusXY, 10*radiusXY, stepZ >
      }  
      finish { diffuse 0.8 ambient 0.2 } 
    }
#end  

union{
  object{ scene_mesh translate <-(maxX)/2.0,-(maxY)/2.0,0 > scale < 1, 1, scaleZ > } 
  // sphere{ <0,0,0>, 0.5*radiusXY }
  texture { tx_plot scale < 1, 1, scaleZ > }
}

#if (closeup)
  // For debugging. Edit these as needed.
  #declare camera_dir = vnormalize(< 00.000, -1.000, +0.500 >);
  #declare camera_disp = 1000*camera_dir;
  #declare camera_radius = 10.0;
  #declare camera_focus = < 0, -0.40*maxY, minZ + 1 >;
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
    rotate -(+50 + viewAzim)*z
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
