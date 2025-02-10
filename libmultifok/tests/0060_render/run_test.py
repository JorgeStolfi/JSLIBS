#! /usr/bin/python3
# Last edited on 2025-02-05 08:35:03 by stolfi

import os, sys, subprocess
from subprocess import run
  
err = sys.stderr

def main():

  run_prog = True

  pixSampling = 2
  
  NX = 120
  NY = 120

  pixDeb = choose_pixDeb(NX, NY)
  outFolder = "out"
  if run_prog:
    clean_output_files(outFolder)
    run_program(NX, NY, pixSampling, outFolder, pixDeb)
  show_results(outFolder,pixDeb)
  os.system(f"( cd out/ && rmempty.sh )");
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
  
  prog = 'test_mfok_render'
  
  err.write("compiling program ...\n")
  status = os.system(f"make {prog}");
  assert status == 0, f"** 'make {prog}' returned status = {status}"

  if not os.path.exists(prog):
    err.write(f"** 'make {prog}' failed\n")
  else:
    err.write(f"running {prog} ...\n")
    status = run(
        [ f"{prog}",
          "-imageSize",   f"{NX}", f"{NY}",
          "-pixSampling", f"{pixSampling}",
          "-debugPixel",  f"{pixDeb[0]}", f"{pixDeb[1]}",
          "-outFolder",   f"{outFolder}"
        ]
      )
    assert status.returncode == 0, f"** '{prog}' failed - returned status = {status}"
  return
  # ------------------------------------------------------------

def show_results(outFolder, pixDeb):  
  # Runs 'display' with all rendered frame images.
  #
  # Assumes that each frame image file is is "{frameFolder}/{imageName}.png"
  # where {frameFolder} is each string returned by {list_frames}.
  
  err.write("collecting the images ...\n")
  images = []
  for frameSub in list_frames(outFolder):
    err.write(f"  frame {frameSub} ...\n");
    images.append(f"{frameSub}/sVal.png")
    images.append(f"{frameSub}/sNrm.png")
    images.append(f"{frameSub}/hAvg.png")

  assert len(images) > 0, "** no images generated"
  # Display the images>
  args = " ".join(images)
  err.write(f"displaying the images from {outFolder} ...\n")
  err.write("args = «" + args + "»\n")
  status = os.system(f"( cd {outFolder} && display -title '%d' -filter Box -resize 'x800<' {args} )")
  return
  # ------------------------------------------------------------
  
def list_frames(outFolder):
  # Returns a list of the names of all frame folders in directory
  # {outFolder}.  The folder names are assumed to begin with "nb".
  
  TS = []
  LS = os.listdir(outFolder); LS.sort()
  FS = [ x for x in LS if x[0:2] == "nb" ]
  for sizeFolder in FS:
    LSS = os.listdir(f"{outFolder}/{sizeFolder}"); LSS.sort()
    FSS = [ f"{sizeFolder}/{x}" for x in LSS if x[0:3] == "lit" ]
    TS += FSS
  return TS
  
main()
