// Last edited on 2012-01-21 12:44:06 by stolfi
    
#macro lamp_array(nrings,skip,angrad,lux) 
  // A round lamp array with {nrings} rings of lamps,
  // (after skipping the first {skip} rings),
  // total angular radius {angrad}, and total strength {lux}.
  // The array is perpendicular to the X axis centered on <1,0,0>.
  // from the origin.
  
  #local totrings = skip+nrings; // Number of rings including skipped ones.
  
  // Compute normalization factor {st_norm} for ring strengths:
  #local k = skip;  // Light ring index.
  #local st_norm = 0;
  #while (k < totrings)
    #local st = lamp_array_rel_ring_strength(k,totrings);
    #local st_norm = st_norm + st; 
    #local k = k + 1;
  #end
  
  #local coin = seed(4615);
  #local totlights = 0;
  union{
    #local k = skip;  // Light ring index.
    #while (k < totrings)
      // Number of lamps {np} in ring {k} and initial phase {phas}:
      #local np = (k = 0 ? 1 : 6*pow(2,k-1));      // Number of lamps in ring.
      #local phas = (k = 0 ? 0 : 0.6180/np);       // Rel. position of first lamp
      // Total relative strength {st} of lamps in ring:
      #local st = lamp_array_rel_ring_strength(k,totrings);
      // Generate the lamps in the ring:
      #local p = 0; // Lamp index.
      #while (p < np)
        // Jittering terms in {k,p}:
        #local dk = (k = 0 ? 0 : 0.20*(2*rand(coin)-1));
        #local dp = (k = 0 ? 0 : 0.20*(2*rand(coin)-1));
        // Angular radius of ring:
        #local tau = lamp_array_rel_ring_radius(k+dk,totrings); // Relative radius of ring
        #local phi = (k = 0 ? 0 : angrad*tau);  // Angular radius of ring.
        #local tht = 360*(phas + (p+dp)/np); 
        light_source {
          <1,0,0>
          color rgb (st/st_norm/np)*lux*<1.0,1.0,1.0>
          rotate phi*z
          rotate tht*x
        }
        #local totlights = totlights + 1;
        // #warning concat("light at phi = ", str(phi,6,1), " tht = ", str(tht,6,1), "\n")
        #local p = p + 1;
      #end
      #local k = k + 1;
    #end
    #debug concat("total lights = ", str(totlights,1,0), "\n")
  }
#end

#macro lamp_array_rel_ring_strength(k,totrings)
  // Relative total strength (unnormalized) of lights in ring {k} of {totrings} rings.
  #local st = (k = 0 ? 0.5*0.5 : 2*k)*((totrings-0.5)*(totrings-0.5) - k*k);
  st
#end

#macro lamp_array_rel_ring_radius(k,totrings)
  // Relative radius of light ring {k} in a lamp cluster with {totrings} rings.
  #local rr = (k = 0 ? 0 : pow(k/(totrings - 0.5),1.5));
  rr
#end
