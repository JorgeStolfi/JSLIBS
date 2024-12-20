#! /usr/bin/python3
# Last edited on 2024-12-17 10:56:09 by stolfi

from math import sin, cos, log, exp, pi, sqrt
import os
import subprocess
import sys

err = sys.stderr

def main():
  light_dirs = (
    ( -0.563, +0.499, +0.658 ),
    ( +0.229, +0.689, +0.686 ),
    ( +0.733, +0.214, +0.645 ),
    ( +0.567, -0.554, +0.608 ),
    ( -0.191, -0.794, +0.576 ),
    ( -0.708, -0.226, +0.668 ),
    #
    ( -0.470, +0.836, +0.282 ),
    ( +0.890, +0.418, +0.185 ),
    ( +0.980, -0.128, +0.151 ),
    ( +0.472, -0.845, +0.251 ),
    ( -0.854, -0.484, +0.190 ),
    ( -0.982, +0.105, +0.155 ),
  )
  run_prog = True      # Shall we run the program?
  scene_type = 'F'     # 'R' is ramp-only, 'F' non-overlapping objs, 'T' overlapping objs.
  single_frame = False # Create a single frame, not a full stack?
  single_light = True  # Consider a single light source, not all light sources?

  ambient = 0.0        # Fraction of light that is non-directional.
  HS = 1               # Use {(2*HS+1)^3} sampoints per pixel.
  KR = 10              # Throw {KR} rays per sampoint. 
  
  scene_WX = 160       # {X}-size of scene in Scene units.
  scene_WY = 120       # {Y}-size of scene in Scene units.
  scene_zMin = 30      # Min {Z} coord of scne in Scene units.
  scene_zMax = 90      # Max {Z} coord of scne in Scene units.
  zDep = 12.0          # Depth of focus in Scene units.
   
  ppu = 2              # Image pixels per Scene unit.
  NX = ppu * scene_WX  # Actual image {X} size.
  NY = ppu * scene_WY  # Actual image {Y} size.
  
  texture = 'noise01'  # See other textures in the "in" folder.
    
  if single_light:
    nLit = 1
  else:
    nLit = len(light_dirs)
    
  patternDir = "in"
  sceneTag = f"st{scene_type}-{texture}"
  sizeTag = f"{NX:04d}x{NY:04d}"
  samplingTag = f"hs{HS:02d}-kr{KR:02d}"
  images = []
  runDir = f"out/{sceneTag}/{sizeTag}-{samplingTag}"
  if run_prog:
    err.write(f"cleaning output folder {runDir} ...\n")
    subprocess.run(["rm", "-rf", f"'{runDir}'"])
  for kLit in range(nLit):
    lid = light_dirs[kLit]
    lightTag = f"L{kLit:03d}-amb{ambient:06.4f}"
    stackDir = f"{runDir}/{lightTag}"

    if single_frame:
      zFoc_min = (scene_zMin + scene_zMax)/2
      zFoc_max = zFoc_min
      zFoc_step = 0.1 # Irrelevant, but cannot be zero.
    else:
      zFoc_step = zDep/2   # Increment of {zFoc} between frames.
      zFoc_min = scene_zMin - zFoc_step/2;
      zFoc_max = scene_zMax + zFoc_step/2;

    if run_prog:
      sys.stderr.write(f"creating frame images for light {lightTag} ...\n")
      subprocess.run(["mkdir", "-pv", f"{stackDir}"])
      prog = "test_mfok_scene_make_frame"
      result = subprocess.run \
        ( [ prog, 
            "-imageSize",    f"{NX}", f"{NY}",
            "-sceneType",    f"{scene_type}",
            "-sceneSize",    "0", f"{scene_WX}", "0", f"{scene_WY}", f"{scene_zMin}", f"{scene_zMax}",
            "-pixSampling",  f"{HS}",
            "-dirSampling",  f"{KR}",
            "-focusHeight",  f"{zFoc_min:8.4f}", "to", f"{zFoc_max:8.4f}", "step", f"{zFoc_step:8.4f}",
            "-dephOfFocus",  f"{zDep}",
            "-patternFile",  f"{patternDir}/{texture}.png",
            "-lightDir",     f"{lid[0]:+8.5f}", f"{lid[1]:+8.5f}", f"{lid[2]:+8.5f}",
            "-ambient",      f"{ambient:6.4f}",
            "-stackDir",     f"{stackDir}",
          ] 
        )
      assert result.returncode == 0, f"** {prog} failed - returned status = {status}"
     
    for fr in list_frames(stackDir) + ("sharp",):
      images.append(f"{lightTag}/{fr}/sVal.png")
    
  err.write("displaying the rendered images ...\n");
  result = subprocess.run \
    ( [ f"( cd {runDir} && display -title '%d/%f' -filter Box -resize '400%' " + \
        " ".join(images) + " )" 
      ], 
      shell=True
    )
  print(result)

  sys.stderr.write("done\n") 
  return 0
    
def list_frames(stackDir):
  # Returns a list of the names of all frame sub-folders in directory
  # {stackDir}, except the "sharp" frame.  The folder names are assumed
  # to begin with "zf-".
  
  LS = os.listdir(stackDir); LS.sort()
  FS = [ x for x in LS if x[0:2] == "zf" ]
  return tuple(FS)

# ----------------------------------------------------------------------
# SHARP_IMAGE := out/${STACKDIR}/frame-sharp/sVal.png
# PIX_PLOT_FILE := out/${STACKDIR}/pixplot.txt
# ----------------------------------------------------------------------
  
main()

