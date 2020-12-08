#! /usr/bin/gawk -f
# Last edited on 2005-12-11 05:29:29 by stolfi

BEGIN {
  set_matrix_YIQ();
  test_cube("YIQ");
  
  set_matrix_YUV();
  test_cube("YUV");
  
  set_matrix_YUV_b();
  test_cube("YUV_b");
  
  set_matrix_YUV_a();
  test_cube("YUV_a");
  
  set_matrix_YCbCr_601_1();
  test_cube("YCbCr_601_1");
  
  fflush("/dev/stdout");
}
          
function set_matrix_YUV(  m)
{
  # European YUV standard
  RY = +0.298911; GY = +0.586611; BY = +0.114478;
  RU = -0.147000; GU = -0.289000; BU = +0.436000;
  RV = +0.615000; GV = -0.515000; BV = -0.100000;
  
  ybias = 0.05;
}
          
function set_matrix_YIQ(  m)
{
  # American YIQ standard

  RY = +0.298911; GY = +0.586611; BY = +0.114478;
  RU = +0.596018; GU = -0.274147; BU = -0.321871;
  RV = +0.211646; GV = -0.522976; BV = +0.311330;
 
  ybias = 0.05;
}
          
function set_matrix_YCbCr_601_1(  m)
{
  m = 65535.0;
  
  RY =  19595.0/m;  GY =   38469.0/m; BY =    7471.0/m;
  RU = -11056.0/m;  GU =  -21712.0/m; BU =   32768.0/m;
  RV =  32768.0/m;  GV =  -27440.0/m; BV =   -5328.0/m;
  
  ybias = 0.05;
}

function set_matrix_YUV_a(  m)
{
  # Original YUV matrix in frgb_to_yuv
  m = 1024;
  
  RY =  306.0/m; GY =  601.0/m; BY =  117.0/m;
  RU =  136.0/m; GU = -136.0/m; BU =    0.0/m;
  RV = -169.0/m; GV =   39.0/m; BV =  130.0/m;
  
  ybias = 0.05;
}
          
function set_matrix_YUV_b(  h)
{
  RY = +0.298911; GY = +0.586611; BY = +0.114478;
  RU = +0.284000; GU = -0.284000; BU = +0.000000;
  RV = -0.139000; GV = -0.085000; BV = +0.224000;
  
  ybias = 0.05;
}

function test_cube(msg,  R,G,B,Y,U,V,d,u,v,h,Um,Vm,UVcos)
{
  printf "----------------------------------------\n";
  printf "%s\n", msg;
  printf "\n";
  
  printf "ybias = %+7.4f\n", ybias;
  Um = sqrt(RU*RU+GU*GU+BU*BU);
  printf "U = [%+7.4f %+7.4f %+7.4f] * %7.4f\n", RU/Um, GU/Um, BU/Um, Um;
  Vm = sqrt(RV*RV+GV*GV+BV*BV);
  printf "V = [%+7.4f %+7.4f %+7.4f] * %7.4f\n", RV/Vm, GV/Vm, BV/Vm, Vm;
  UVcos = (RU*RV + GU*GV + BU*BV)/(Um*Vm);
  printf "cos(U*V) = %7.4f\n", UVcos;

  # Get max chroma (assume it comes from RGB=(0,0,1)):
  h = sqrt(BU*BU + BV*BV)*(1.0 + ybias)/(BY + ybias);
  printf "h = %9.6f\n", h;

  printf "\n";
  
  for (R = 0; R <= 1; R++)
    for (G = 0; G <= 1; G++)
      for (B = 0; B <= 1; B++)
        { 
          Y = R*RY + G*GY + B*BY;
          U = R*RU + G*GU + B*BU;
          V = R*RV + G*GV + B*BV;
          
          d = (1.0 + ybias)/(Y + ybias)/h;
          u = U*d;
          v = V*d;
          
          printf " RGB=[%+6.3f %+6.3f %+6.3f]", R, G, B; 

          printf " YUV=[%+6.3f %+6.3f %+6.3f]", Y, U, V; 
          printf " C=%5.3f", sqrt(U*U+V*V); 
          
          printf " d=%7.3f", d; 
          
          printf " Yuv=[%+6.3f %+6.3f %+6.3f]", Y, u, v; 
          printf " c=%5.3f", sqrt(u*u+v*v); 

          printf "\n";
        }
}
