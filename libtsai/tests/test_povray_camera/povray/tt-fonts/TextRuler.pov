/*
   Text Ruler  V.02021997
   A tool to determine the length of a text string in POV-Ray
   Written by Patrick Bass

   -----------------------------------------------------------

   Type in the string you need to measure into TEXT.
   Set the path and font name to FONT.
   Render.
   Read the string length right off the ruler.

   Accurate to about 1/80 of a unit.
*/




// STEP 1:  ENTER YOUR TEXT ON THE LINE BELOW.

#declare TEXT="POV-Ray"

// STEP 2:  ENTER WHERE AND WHAT YOUR FONT IS.
// Typically like "C:\WINDOWS\SYSTEM\WINGDINGS.TTF"

#declare FONT="\winnt\fonts\arial.ttf"

// STEP3:  RENDER AN OVERVIEW
// Declare VIEW=1 to render an overview so you can
// get an idea how big the text is.

// Finally, when you have an idea of about where the text ends
// from the overview (like: about 4 something)
// set CAM_X to equal that number, then
// Declare VIEW=2 and render.
// Read the length off the ruler.

// EACH SQUARE ON THE FLOOR IS ONE POV UNIT SQUARE.
// So the word "POV-Ray" looks to be about four units
// wide in the overview.

#declare VIEW=1

#declare CAM_X=4.15








 //  User modification below not really needed,
 //  but play with it to your hearts content!
 //  --------------------------------------------------


#include "colors.inc"
#include "woods.inc"
#include "metals.inc"
#include "golds.inc"

text{
   ttf FONT TEXT, .1, 0
   texture { pigment { ForestGreen } finish { reflection 0.1 }  }
}

#if (VIEW=1)
   #declare CAMERA_LOCATION=< 6, 8, -8>
   #declare LOOK_AT=<6, 0, 0>
#else
   #declare CAMERA_LOCATION=< CAM_X, 0.5, -1>
   #declare LOOK_AT=< CAM_X, 0.25, 0>
#end

camera {
   location CAMERA_LOCATION
   look_at LOOK_AT
}

light_source { <20,20,-20> color White*2 }

#declare HRULERHEIGHT=0.025
#declare HRULERWIDTH = (-0.25)

#declare HRULERWOOD=
object {
   difference {
      box {
         <0, 0, 0>,<12, HRULERHEIGHT, HRULERWIDTH>
         texture { T_Wood11 scale 0.5 rotate 90*y }
      }
      cylinder {
         <-1, 0, 0>, <13, 0, 0>, 0.02
         texture { T_Wood11 scale 0.5 rotate 90*y }
         translate -0.125*z
         translate (HRULERHEIGHT+0.01)*y
      }
   }
}


#declare TICK=
object {
   box { <0, 0, 0>,<0.005, HRULERHEIGHT+0.005, -0.075>
      texture { T_Gold_3C }
   }
}

#declare INCHES=(0)
#declare XP=(0)


#declare HTICKMARKS=
object {
   union {
      #while (INCHES<12)
         #declare XP=(0)
         #while (XP<10)

            #declare TICKPLACE= (x*(INCHES+(XP/10)+0.0001) )
            object { TICK
               translate TICKPLACE
            }
            object { TICK
               scale <1, 1, 0.5>
               translate x*(TICKPLACE+0.05)
            }
            object { TICK
               scale <1, 1, 0.25>
               translate x*(TICKPLACE+0.025)
            }
            object { TICK
               scale <1, 1, 0.25>
               translate x*(TICKPLACE+0.075)
            }

            #if (INCHES>0)
               #if (INCHES>9)
                  object {
                     text{
                        ttf FONT
                        concat ("1",chr(48+(INCHES-10))), .1, 0
                        scale 0.05
                        pigment { Black }
                     }
                     rotate 90*x
                     translate x*(INCHES+(XP/10))
                     translate (HRULERHEIGHT*y)+0.001
                     translate z*(-0.2)
                  }
               #else
                  object {
                     text{
                        ttf FONT
                        chr(48+INCHES), .1, 0
                        scale 0.05
                        pigment { Black }
                     }
                     rotate 90*x
                     translate x*(INCHES+(XP/10))
                     translate (HRULERHEIGHT*y)+0.001
                     translate z*(-0.2)
                  }
               #end
            #end
            object {
               text{
                  ttf FONT
                  chr(48+XP), .1, 0
                  scale 0.05
                  pigment { Black }
               }
               rotate 90*x
               translate x*(INCHES+(XP/10))
               translate (HRULERHEIGHT*y)+0.001
               translate z*(-0.25)
            }
            #declare XP=(XP+1)
         #end
         #if (INCHES=1)
            object {
               text{
                  ttf FONT concat(chr(48+INCHES)," unit"), .1, 0
                  scale 0.15
                  texture { T_Gold_3C }
               }
         #else
            object {
               text{
                  ttf FONT concat(chr(48+INCHES)," units"), .1, 0
                  scale 0.15
                  texture { T_Gold_3C }
               }
         #end
               rotate 90*x
               translate x*INCHES
               translate .01*y
               translate z*(-0.4)
            }
            #declare INCHES=(INCHES+1)
      #end
   }
}



plane { y, 0
   pigment {
      checker color White color Gray90
      scale 1
   }
}
plane { z, 0.05
   pigment {
      checker color White color Gray90 scale .05
   }
}

#declare HRULER=
object {
   union {
      object { HRULERWOOD }
      object { HTICKMARKS }
   }
}

object { HRULER }


//end of source code

