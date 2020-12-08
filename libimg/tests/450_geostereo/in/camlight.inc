// Câmeras e luzes padronizadas
// Last edited on 2017-04-24 15:57:31 by stolfilocal

#macro camlight(ctr,rad,cav,dst,upp,lux)
  // {ctr} = center of interest in scene.
  // {rad} = approx radius of scene.
  // {cav} = vector pointing from {ctr} to camera; only its direction matters.
  // {dst} = the distance from camera to {ctr}.
  // {upp} = the scene's vertical axis ({sky} parameter), usually {z} or {y}.
  // {lux} = a scaling factor for the intensity of standard light sources.
  #local dir = vnormalize(cav);
  #local len = (dst);
  #local swh = sqrt(image_width/image_height);
  #local ape = 1.41*(rad)/(len);  // Camera aperture.
  camera {
    location (ctr) + len*dir
    right  -swh*ape*x
    up     1.0/swh*ape*y
    sky    (upp)
    look_at (ctr)
  }
  #if (lux > 0)
    #local cam_elev = degrees(atan2(dir.z,sqrt(dir.x*dir.x + dir.y*dir.y)));
    #local cam_azim = degrees(atan2(dir.y,dir.x));
    #local nrings = 2;
    #local angrad = 25;
    object{ lamp_array(nrings,angrad,lux)
      rotate 360*(sqrt(5)-1)/2*x
      scale 10*len
      rotate -(90 - 0.667*(90 - cam_elev))*y 
      rotate (cam_azim + 15)*z
      translate ctr 
    }
  #end
#end
    
#macro lamp_array(nrings,angrad,lux) 
  // A round lamp array with {nrings} rings of lamps,
  // total angular radius {angrad}, and total strength {lux}.
  // The array is centered on the X axis at distance 1 
  // from the origin.

  // Compute normalization factor {st_norm} for ring strengths:
  #local k = 0;  // Light ring index.
  #local st_norm = 0;
  #while (k < nrings)
    #local st = camlight_rel_ring_strength(k,nrings);
    #local st_norm = st_norm + st; 
    #local k = k + 1;
  #end
  
  #local coin = seed(4615);
  
  union{
    #local k = 0;  // Light ring index.
    #while (k < nrings)
      // Number of lamps {np} in ring {k} and initial phase {phas}:
      #local np = (k = 0 ? 1 : 6*pow(2,k-1));      // Number of lamps in ring.
      #local phas = (k = 0 ? 0 : 0.6180/np);       // Rel. position of first lamp
      // Total relative strength {st} of lamps in ring:
      #local st = camlight_rel_ring_strength(k,nrings);
      // Generate the lamps in the ring:
      #local p = 0; // Lamp index.
      #while (p < np)
        // Jittering terms in {k,p}:
        #local dk = (k = 0 ? 0 : 0.20*(2*rand(coin)-1));
        #local dp = (k = 0 ? 0 : 0.20*(2*rand(coin)-1));
        // Angular radius of ring:
        #local tau = camlight_rel_ring_radius(k+dk,nrings); // Relative radius of ring
        #local phi = (k = 0 ? 0 : angrad*tau);  // Angular radius of ring.
        #local tht = 360*(phas + (p+dp)/np); 
        light_source {
          <1,0,0>
          color rgb (st/st_norm/np)*lux*<1.0,1.0,1.0>
          rotate phi*z
          rotate tht*x
        }
        // #warning concat("light at phi = ", str(phi,6,1), " tht = ", str(tht,6,1), "\n")
        #local p = p + 1;
      #end
      #local k = k + 1;
    #end
  }
#end

#macro camlight_rel_ring_strength(k,nrings)
  // Relative total strength (unnormalized) of lights in ring {k} of {nrings} rings.
  #local st = (k = 0 ? 0.5*0.5 : 2*k)*((nrings-0.5)*(nrings-0.5) - k*k);
  st
#end

#macro camlight_rel_ring_radius(k,nrings)
  // Relative radius of light ring {k} in a lamp cluster with {nrings} rings.
  #local rr = (k = 0 ? 0 : pow(k/(nrings - 0.5),1.5));
  rr
#end
