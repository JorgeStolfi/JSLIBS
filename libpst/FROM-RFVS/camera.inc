//  Camera placement for the MP logo
//  Last edited on 2007-01-05 23:14:24 by stolfi

//  Requires 
//    camera_focus  = point of interest in scene.
//    camera_radius = nominal radius of scene around {camera_focus}.
//    camera_disp   = displacement vector from {camera_focus} to camera position.
//    image_rel_{width,height} = relative dimensions of image.

//  Defines:
//    The camera.
//    camera_pos = position of camera.
//    camera_distance = distance from camera to focus.
//    camera_dir = unit vector from focus towards camera.
//    camera_{h,v}_dir = unit world vectors parallel to sides of image (right and up).
//    camera_{h,v}_{min,max} = horizontal and vertical extent of image at focus.

#local camera_pos = camera_focus + camera_disp;
#declare camera_distance = vlength(camera_disp);
#declare camera_dir = camera_disp/camera_distance;

#local sky_dir = z;
#declare camera_h_dir = vnormalize(vcross(sky_dir, camera_dir));
#declare camera_v_dir = vnormalize(vcross(camera_dir, camera_h_dir));

#local slop = 0.77;

#local c_h = sqrt(image_rel_width/image_rel_height);
#local c_v = sqrt(image_rel_height/image_rel_width);

#declare camera_h_min = - slop * c_h * camera_radius;
#declare camera_h_max = + slop * c_h * camera_radius;
#declare camera_v_min = - slop * c_v * camera_radius;
#declare camera_v_max = + slop * c_v * camera_radius;

#local c_a_h = (camera_h_max - camera_h_min)/camera_distance;
#local c_a_v = (camera_v_max - camera_v_min)/camera_distance;

camera {
  location  camera_pos
  right     c_a_h * <00, -1, 00>
  up        c_a_v * <00, 00, +1>
  sky       sky_dir
  look_at   camera_focus
} 
