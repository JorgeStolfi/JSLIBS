// Câmeras e luzes padronizadas
// Last edited on 2012-01-21 12:41:16 by stolfi

#macro camlight(ctr,rad,cav,dst,upp,lux,lnr,lar,lct)
  // {ctr} = center of interest in scene.
  // {rad} = approx radius of scene.
  // {cav} = vector pointing from {ctr} to camera; only its direction matters.
  // {dst} = the distance from camera to {ctr}.
  // {upp} = the scene's vertical axis ({sky} parameter), usually {z} or {y}.
  // {lux} = a scaling factor for the intensity of standard light sources.
  // {lnr} = number of rings in light cluster.
  // {lar} = angular radius of light cluster.
  // {lct} = 0 place light away from camera, 1 place around camera (flash-like).
  camera_standard(ctr,rad,cav,dst,upp)
  light_default(ctr,rad,cav,dst,upp,lux,lnr,lar,lct)
#end

#macro camera_standard(ctr,rad,cav,dst,upp)
  // {ctr} = center of interest in scene.
  // {rad} = approx radius of scene.
  // {cav} = vector pointing from {ctr} to camera; only its direction matters.
  // {dst} = the distance from camera to {ctr}.
  // {upp} = the scene's vertical axis ({sky} parameter), usually {z} or {y}.
  #local dir = vnormalize(cav);
  #local swh = sqrt(image_width/image_height);
  #local ape = 1.41*(rad)/(dst);  // Camera aperture.
  camera {
    location (ctr) + (dst)*dir
    right  -swh*ape*x
    up     1.0/swh*ape*y
    sky    (upp)
    look_at (ctr)
  }
#end

#macro light_default(ctr,rad,cav,dst,upp,lux,lnr,lar,lct)
  // {ctr} = center of interest in scene.
  // {rad} = approx radius of scene.
  // {cav} = vector pointing from {ctr} to camera; only its direction matters.
  // {dst} = the distance from camera to {ctr}.
  // {upp} = the scene's vertical axis ({sky} parameter), usually {z} or {y}.
  // {lux} = a scaling factor for the intensity of standard light sources.
  // {lnr} = number of rings in light cluster.
  // {lar} = angular radius of light cluster.
  // {lct} = 0 place light away from camera, 1 place around camera (flash-like).
  
  // !!! Should account for {upp}.
  #local dir = vnormalize(cav);
  #if (lux > 0)
    #local dir_xy = vlength( < dir.x, dir.y, 0 > );
    #if (abs(dir.z) > 100*dir_xy)
      // Camera is nearly vertical.
      #local cam_azim = 0;
      #if (dir.z > 0)
        #local cam_elev = 90;
      #else
        #local cam_elev = -90;
      #end
    #else
      // Camera is not very vertical.
      #local cam_elev = degrees(atan2(dir.z,dir_xy));
      #local cam_azim = degrees(atan2(dir.y,dir.x));
    #end
    
    #if (lct > 0)
      #local lux_elev = cam_elev;
      #local lux_azim = cam_azim;
    #else
      #local lux_elev = cam_elev + 15;
      #local lux_azim = cam_azim + 15;
    #end
  
    object{ lamp_array(lnr,1,lar,lux)
      scale 10*(dst)
      rotate 360*(sqrt(5)-1)/2*x
      rotate -lux_elev*y
      rotate lux_azim*z
      translate ctr 
    }
  #end
#end
