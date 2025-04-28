#! /usr/bin/python3
# Last edited on 2025-04-19 06:45:24 by stolfi

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
  singleLight = False # Use a single light, not all light?

  if singleLight:
    iniLight = 0                # Index of first light source to use
    finLight = 0                # Index of last light source to use.
  else:
    iniLight =  0               # Index of first light source to use.
    finLight = len(lightDirs)-1 # Index of last light source to use.

  HS = 1              # Use {(2*HS+1)^2} sampoints per pixel.
  KR = 7              # Throw {KR} rays per sampoint. 
  
  images = [];
  for lightNum in range(iniLight, finLight+1):
    images += make_stack \
      ( sceneType = 'Q', 
        lightNum = lightNum, lightDir = lightDirs[lightNum],
        ambient = 0.00, gloss = 0.00,
        runProg = runProg, singleFrame = singleFrame,
        HS = HS, KR = KR
      )
      
  show_images(images)
  return
  # ......................................................................
  
def make_stack(sceneType, lightNum, lightDir, ambient, gloss, runProg, singleFrame, HS, KR):  
  # Runs the test program with specified parameters.
  # 
  #   {sceneType} is a letter 'R', 'S' etc (see below)
  #   {lightNum} is the light source index, for file names (see below)
  #   {lightDir} is a 3-list, the coords of the unit vector pointing towards the light source.
  #   {ambient} is the relative amount of ambient (isotropic) light.
  #   {gloss} is the relative amount of glossy component in the BRDF.
  #   {runProg} is true to run the program, false to only display results of a previous run.
  #   {HS} pixel subsampling param; will use {(2*HS+1)^2} sampoints per pixel.
  #   {KR} ray sampling param; will throw {KR} rays per sampoint. 

  # The procedure writes a stack of frames to a folder {stackFolder}
  # that is "out/{sceneTag}/{sizeTag}-{samplingTag}/{lightTag}" where
  
  #   {sceneTag} is "st{sceneType}-{texture}-amb{A.AA}-glo{G.GG}"
  #   {sizeTag} is "{NX}x{NY}", both formatted as '%04d'
  #   {samplingTag} is "hs{HH}-kr{RR}"
  #   {lightTag} is "L{LLLL}"
  #   {texture} is the texture file, e.g. "noise01"
  #   {HH} is the pixel subsampling parameter {HS} formatted as '%02d' 
  #   {RR} is the ray per sampoint parameter {KR} as '%02d' 
  #   {LLLL} is the light source index as '%03d'
  #   {A.AA} is the relative amount of isotropic ("ambient") light as '%4.2f'
  #   {G.GG} is the relative amount of glossy surface finish as '%4.2f'

  # Scene types:
  #  'D' Single disk.
  #  'B' Single ball.
  #  'C' Single cone.
  #  'P' Single square pyramid.
  
  #  'R' Ramp only.
  #  'F' Non-overlapping mixed objects.
  #  'T' Overlapping mixed objects.
  #  'Q' Four non-overlapping objects, one of each type.
  
  # SCU = Scene unit of length.
  scene_WX = 512       # {X}-size of scene in SCU.
  scene_WY = 384       # {Y}-size of scene in SCU.
  scene_zMin = 20      # Min {Z} coord of scene in SCU.
  scene_WXY = min(scene_WX, scene_WY);   # Usefuel {XY] size of scene in SCU. 
  # Define max scene {Z} coord {scene_zMax} in SCU:
  if sceneType == 'R' or sceneType == 'F' or sceneType == 'T':
    scene_zMax = scene_zMin + 50
  elif sceneType == 'Q':
    scene_zMax = scene_zMin + 2*scene_WXY/3
  else:
    scene_zMax = scene_zMin + scene_WXY
    
  NF = 12              # Num of frames in full stack.
 
  zFoc_step = (scene_zMax - scene_zMin)/(NF-2)  # Inc of {zFoc} btw frames in SCU.  
  zDep = 2*zFoc_step   # Depth of focus in SCU.
   
  ppu = 1              # Image pixels per Scene unit.
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
  Er("lightTag = %s\n" % lightTag)
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
    Er(f"creating multi-focus stack for light {lightTag} ...\n")
    Er(f"cleaning output folder {stackFolder} ...\n")
    bash(f"rm -rf {stackFolder}")
    bash(f"mkdir -pv {stackFolder}")
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
    Er("command = [ " + " ".join(command) + " ]\n");
    run_command(command);

  for fr in list_frames(stackFolder) + ("sharp",):
    images.append(f"{stackFolder}/{fr}/sVal.png")

  Er("done\n") 
  return images
  
def show_images(images):
  Er("displaying the rendered images ...\n");
  bash(f"( display -title '%d/%f' -filter Box -resize '200%' " + " ".join(images) + " )")
  return
  # ......................................................................
    
def list_frames(stackFolder):
  # Returns a list of the names of all frame sub-folders in directory
  # {stackFolder}, except the "sharp" frame.  The folder names are assumed
  # to begin with "zf-".
  
  LS = os.listdir(stackFolder); LS.sort()
  FS = [ x for x in LS if x[0:2] == "zf" ]
  return tuple(FS)
  # ......................................................................

def file_OK(path):
  return path != None and os.path.exists(path) and os.path.getsize(path) > 0
  # ......................................................................

def Er(msg):
  sys.stderr.write(msg)
  # ......................................................................

def run_command(command):
  result = subprocess.run(command, text = True)
  if result.returncode != 0:
    print(result.stderr)
    print(result.stdout)
    assert False, f"** {command[0]} failed - returned status = {result}"
  return
  # ......................................................................
  
def bash(cmd):
  # Execute the string {cmd} with "/bin/bash".
  result = subprocess.run([ cmd ], shell = True, executable = "/bin/bash")
  if result.returncode != 0:
    print(result.stderr)
    print(result.stdout)
    assert False, f"** {cmd} failed - returned status = {result}"
  # ......................................................................

# ----------------------------------------------------------------------
# SHARP_IMAGE := out/${STACKFOLDER}/frame-sharp/sVal.png
# PIX_PLOT_FILE := out/${STACKFOLDER}/pixplot.txt
# ----------------------------------------------------------------------
  
main()

