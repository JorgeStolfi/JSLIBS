#! /usr/bin/python3
# Last edited on 2025-03-08 15:38:36 by stolfi

from math import sin, cos, log, exp, pi, sqrt
import os
import subprocess
import sys

err = sys.stderr

def main():
  # This script calls {test_mfok_scene_make_frame} for various
  # lighting and finish parameters.

  lightDirs = (
    ( +0.567, -0.554, +0.608 ),
    ( -0.191, -0.794, +0.576 ),
    ( -0.708, -0.226, +0.668 ),
    ( -0.563, +0.499, +0.658 ),
    ( +0.229, +0.689, +0.686 ),
    ( +0.733, +0.214, +0.645 ),
    #
    ( +0.472, -0.845, +0.251 ),
    ( -0.854, -0.484, +0.190 ),
    ( -0.982, +0.105, +0.155 ),
    ( -0.470, +0.836, +0.282 ),
    ( +0.890, +0.418, +0.185 ),
    ( +0.980, -0.128, +0.151 ),
  )

  runProg = True      # Actually run the program?
  singleFrame = False # Create a single frame, not a full stack?

  iniLight =  0       # Index of first light source to use.
  finLight = 11       # Index of last light source to use.
  
  images = [];
  for lightNum in range(iniLight, finLight+1):
    images += make_stack \
      ( sceneType = 'Q', 
        lightNum = lightNum, lightDir = lightDirs[lightNum],
        ambient = 0.00, gloss = 0.00,
        runProg = runProg, singleFrame = singleFrame
      )
      
  show_images(images)
  
def make_stack(sceneType, lightNum, lightDir, ambient, gloss, runProg, singleFrame):  
  # Runs the program, that writes a stack 
  # of frames to a folder {stackFolder} that is
  # "out/{sceneTag}/{sizeTag}-{samplingTag}/{lightTag}"
  # where 
  
  #   {sceneTag} is "st{sceneType}-{texture}-amb{A.AA}-glo{G.GG}"
  #   {sizeTag} is "{NX}x{NY}", both formatted as '%04d'
  #   {samplingTag} is "hs{HH}-kr{RR}"
  #   {lightTag} is "L{LLLL}"
  # 
  #   {sceneType} is a letter 'R', 'S' etc (see below)
  #   {texture} is the texture file, e.g. "noise01"
  #   {HH} is the pixel subsampling parameter {HS} formatted as '%02d' 
  #   {RR} is the number of rays per subsampling point formatted as '%02d' 
  #   {LLLL} is the light source index formatted as '%03d'
  #   {A.AA} is the relative amount of isotropic ("ambient") light formatted as '%4.2f'
  #   {G.GG} is the relative amount of glossy surface finish formatted as '%4.2f'

  # Scene types:
  #  'D' Single disk.
  #  'B' Single ball.
  #  'C' Single cone.
  #  'P' Single square pyramid.
  
  #  'R' Ramp only.
  #  'F' Non-overlapping mixed objects.
  #  'T' Overlapping mixed objects.
  #  'Q' Four non-overlapping objects, one of each type.

  HS = 1               # Use {(2*HS+1)^2} sampoints per pixel.
  KR = 10              # Throw {KR} rays per sampoint. 
  
  # SCU = Scene unit of length.
  scene_WX = 256       # {X}-size of scene in SCU.
  scene_WY = 192       # {Y}-size of scene in SCU.
  scene_zMin = 20      # Min {Z} coord of scne in SCU.
  scene_WXY = min(scene_WX, scene_WY);   # Usefuel {XY] size of scene in SCU. 
  # Define max scene {Z} coord {scene_zMax} in SCU:
  if sceneType == 'R' or sceneType == 'F' or sceneType == 'T':
    scene_zMax = scene_zMin + 50
  elif sceneType == 'Q':
    scene_zMax = scene_zMin + 100
  else:
    scene_zMax = scene_zMin + scene_WXY
    
  NF = 12              # Num of frames in full stack.
 
  zFoc_step = (scene_zMax - scene_zMin)/(NF-2)  # Inc of {zFoc} btw frames in SCU.  
  zDep = 2*zFoc_step   # Depth of focus in SCU.
   
  ppu = 2              # Image pixels per Scene unit.
  NX = ppu * scene_WX  # Actual image {X} size.
  NY = ppu * scene_WY  # Actual image {Y} size.

  texture = 'melon14'  # See other textures in the "in" folder.
    
  patternFolder = "in"
  sceneTag = f"st{sceneType}-{texture}-amb{ambient:04.2f}-glo{gloss:04.2f}"
  sizeTag = f"{NX:04d}x{NY:04d}"
  samplingTag = f"hs{HS:02d}-kr{KR:02d}"
  images = []
  runFolder = f"out/{sceneTag}/{sizeTag}-{samplingTag}"

  lid = lightDir
  lightTag = f"L{lightNum:03d}"
  err.write("lightTag = %s\n" % lightTag)
  stackFolder = f"{runFolder}/{lightTag}"

  if singleFrame:
    zFoc_min = (scene_zMin + scene_zMax)/2
    zFoc_max = zFoc_min
    # The {zFoc_step} is now irrelevant but {zDep} is not.
  else:
    zFoc_min = scene_zMin - zFoc_step/2;
    zFoc_max = scene_zMax + zFoc_step/2;
    # With {zFoc_step} as defined above, should give {NF} frames.

  if runProg:
    err.write(f"creating multi-focus stack for light {lightTag} ...\n")
    err.write(f"cleaning output folder {stackFolder} ...\n")
    subprocess.run(["rm", "-rf", f"{stackFolder}"])
    subprocess.run(["mkdir", "-pv", f"{stackFolder}"])
    prog = "test_mfok_scene_make_frame"
    command = \
      [ prog, 
        "-imageSize",    f"{NX}", f"{NY}",
        "-sceneType",    f"{sceneType}", 
        "-sceneSize",    "0", f"{scene_WX}", "0", f"{scene_WY}", f"{scene_zMin}", f"{scene_zMax}",
        "-pixSampling",  f"{HS}",
        "-dirSampling",  f"{KR}",
        "-focusHeight",  f"{zFoc_min:8.4f}", "to", f"{zFoc_max:8.4f}", "step", f"{zFoc_step:8.4f}",
        "-dephOfFocus",  f"{zDep}",
        "-patternFile",  f"{patternFolder}/{texture}.png",
        "-lightDir",     f"{lid[0]:+8.5f}", f"{lid[1]:+8.5f}", f"{lid[2]:+8.5f}",
        "-ambient",      f"{ambient:4.2f}",
        "-gloss",        f"{gloss:4.2f}",
        "-stackFolder",  f"{stackFolder}",
      ] 
    err.write("command = [ " + " ".join(command) + " ]\n");
    result = subprocess.run( command );
    assert result.returncode == 0, f"** {prog} failed - returned status = {result}"

  for fr in list_frames(stackFolder) + ("sharp",):
    images.append(f"{stackFolder}/{fr}/sVal.png")

  err.write("done\n") 
  return images
  
def show_images(images):
  err.write("displaying the rendered images ...\n");
  result = subprocess.run \
    ( [ f"( display -title '%d/%f' -filter Box -resize '200%' " + \
        " ".join(images) + " )" 
      ], 
      shell=True
    )
  print(result)
    
def list_frames(stackFolder):
  # Returns a list of the names of all frame sub-folders in directory
  # {stackFolder}, except the "sharp" frame.  The folder names are assumed
  # to begin with "zf-".
  
  LS = os.listdir(stackFolder); LS.sort()
  FS = [ x for x in LS if x[0:2] == "zf" ]
  return tuple(FS)

# ----------------------------------------------------------------------
# SHARP_IMAGE := out/${STACKFOLDER}/frame-sharp/sVal.png
# PIX_PLOT_FILE := out/${STACKFOLDER}/pixplot.txt
# ----------------------------------------------------------------------
  
main()

