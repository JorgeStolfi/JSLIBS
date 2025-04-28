#! /usr/bin/python3
# Last edited on 2025-04-19 06:45:59 by stolfi

import os, sys, subprocess
from subprocess import run
  
err = sys.stderr
  
prog = 'test_mfok_render'

def main():

  runProg = False

  pixSampling = 2
  
  NX = 120
  NY = 120

  pixDeb = choose_pixDeb(NX, NY)
  outFolder = "out"
  if runProg:
    clean_output_files(outFolder)
    run_program(NX, NY, pixSampling, outFolder, pixDeb)
  show_results(outFolder,pixDeb)
  bash(f"( cd out/ && rmempty.sh )");
  return 0
  
def clean_output_files(outFolder):
  # Removes all output files.
  os.system(f"cd {outFolder}/ && rm -rf nb*")
  return
  # ------------------------------------------------------------

def choose_pixDeb(NX, NY):
  # Returns a list {(ix,iy)} with the indices of the pixel 
  # for which the ray-tracing, sampling and debugging data is to 
  # be writen.
  
  # pixDeb = ( NX//2, NY//2 )
  pixDeb = ( 4, 5 )
  return pixDeb
  # ------------------------------------------------------------

def run_program(NX, NY, pixSampling, outFolder, pixDeb):
  # Runs the C program generating the rendered frame images.
  # Also compiles the program.
  
  Er("compiling program ...\n")
  bash(f"make {prog}");

  if not os.path.exists(prog):
    Er(f"** 'make {prog}' failed\n")
  else:
    Er(f"running {prog} ...\n")
    cmd = \
        [ f"{prog}",
          "-imageSize",   f"{NX}", f"{NY}",
          "-pixSampling", f"{pixSampling}",
          "-debugPixel",  f"{pixDeb[0]}", f"{pixDeb[1]}",
          "-outFolder",   f"{outFolder}"
        ]
    run_command(cmd)
  return
  # ------------------------------------------------------------

def show_results(outFolder, pixDeb):  
  # Runs 'display' with all rendered frame images.
  #
  # Assumes that each frame image file is is "{frameFolder}/{imageName}.png"
  # where {frameFolder} is each string returned by {list_frames}.
  
  Er("collecting the images ...\n")
  images = []
  for frameSub in list_frames(outFolder):
    Er(f"  frame {frameSub} ...\n");
    for comp in ( "sVal", "sNrm", "hAvg", ):
      fName = f"{frameSub}/{comp}.png"
      if file_OK(f"{outFolder}/{fName}"):
        images.append(fName)

  assert len(images) > 0, "** no images generated"
  # Display the images>
  args = " ".join(images)
  Er(f"displaying the images from {outFolder} ...\n")
  Er("args = «" + args + "»\n")
  bash(f"( cd {outFolder} && display -title '%d' -filter Box -resize 'x800<' {args} )")
  return
  # ------------------------------------------------------------
  
def list_frames(outFolder):
  # Returns a list of the names of all frame folders in directory
  # {outFolder}, sans the "{outFolder}/" part.  The folder names
  # are assumed to begin with "nb*/lit*".
  
  frameFolders = []
  sizeFolders = os.listdir(outFolder); sizeFolders.sort()
  sizeFolders = [ x for x in sizeFolders if x[0:2] == "nb" ]
  for sizeFolder in sizeFolders:
    litFolders = os.listdir(f"{outFolder}/{sizeFolder}"); litFolders.sort()
    litFolders = [ x for x in litFolders if x[0:3] == "lit" ]
    for litFolder in litFolders:
      frameFolder = f"{sizeFolder}/{litFolder}"
      frameFolders.append(frameFolder)
  return frameFolders

def Er(msg):
  sys.stderr.write(msg)
  # ......................................................................

def file_OK(path):
  return path != None and os.path.exists(path) and os.path.getsize(path) > 0
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

main()
