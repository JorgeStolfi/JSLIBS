global_settings { assumed_gamma 2.2 }

#include "colors.inc"
#include "shapes.inc"

background { color rgb <1.000, 1.000, 1.000>  }
#include "scene.inc"

#declare maxSize = max(maxX,maxY);

object{ 
	union{
		object{ 
		scene_mesh
			translate <-(maxX)/2.0,-(maxY)/2.0,0 >
  			rotate 60*z 
		}
		
	}
	
}

light_source { <20000, 50000, 100000> color rgb 0.80*<1,1,1> }
light_source { <-20000, 50000, 100000> color rgb 0.30*<1,1,1> }

#declare camera_disp = 5*maxSize*<0.3,0.4,0.4>;
#declare camera_radius = 0.7*maxSize;
#declare camera_focus = <0,0,0> ;
#declare image_rel_width = image_width;
#declare image_rel_height = image_height;


#include "camera.inc"