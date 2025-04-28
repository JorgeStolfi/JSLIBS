#! /usr/bin/python3 
# Last edited on 2025-04-27 21:56:05 by stolfi

from math import sin, cos, log, exp, pi, sqrt
import os, sys
import subprocess

PROG = "test_spots_image" 
TEXTURE_BANK = "projects/image-collections/texture-bank"

def main():

  relSpotDist = 4.0/3.0
  maxSpots = 20000
  imageSize = 512
  
  # spotSizes = []
  # for spotSize_max in range(10):
  #   for spotSize_min in range(spotSize_max):
  #     spotSizess.append((spotSize_min, spotSize_max))
  # spotSizes.reverse()
  
  spotSizes = []
  spotSize_var = 3
  for spotSize_min in range(10 - spotSize_var):
    spotSize_max = spotSize_min + spotSize_var
    spotSizes.append((spotSize_min, spotSize_max))
  spotSizes.reverse()

  # spotSizes = ( (0,5), (1,6), (2,7), (3,8), (4,9), )
  # spotSizes = ( (6,9), (5,8), (4,5), )
  # spotSizes = ( (0,5), )
  # spotSizes = ( (4,9), )
  # spotSizes = ( (0,0), )
  
  for spotSize_min, spotSize_max in spotSizes:
    numPasses =  50 + 30*spotSize_max
    run_single(imageSize, imageSize, spotSize_min, spotSize_max, relSpotDist, maxSpots, numPasses)
  
  return
  # ......................................................................

def run_single(NX, NY, spotSize_min, spotSize_max, relSpotDist, maxSpots, numPasses):
  # Runs the program {PROG}.
  outDir = f"out/png-{NX:04d}x{NY:04d}"
  bash(f"mkdir -pv {outDir}")

  imageName = f"spots-{spotSize_min}{spotSize_max}"
  outPrefix = f"{outDir}/{imageName}"

  outBegTotalEnergyPlot = f"{outPrefix}-eplot-beg.txt"
  outEndTotalEnergyPlot = f"{outPrefix}-eplot-end.txt"

  outBegEnergyImage = f"{outPrefix}-espot-beg.png"
  outEndEnergyImage = f"{outPrefix}-espot-end.png"

  outBegSpotsImage = f"{outPrefix}-00000.png"
  outEndSpotsImage = f"{outPrefix}-final.png"

  Er(f"image size = {NX}x{NY} spot sizes = {spotSize_min} .. {spotSize_max}\n")
  Er(f"relSpotDist = {relSpotDist:.3f}\n")
  
  single_clean(outPrefix)
  
  Er(f"=== making {outEndSpotsImage} ...\n")
  command = \
    [ PROG, 
      "-imageSize",    f"{NX}", f"{NY}",
      "-spotSizes",    f"{spotSize_min}",  f"{spotSize_max}", 
      "-relSpotDist",  f"{relSpotDist:.4f}",
      "-maxSpots",     f"{maxSpots:d}",
      "-numPasses",    f"{numPasses:d}",
    ] 
  Er("command = [ " + " ".join(command) + " ]\n");
  run_command( command );

  bash(f"ls -l {outBegSpotsImage} {outEndSpotsImage}");
  
  resize = 200 if NY >= 200 else 400 # Percentage resize factor for display.

  Er(f"=== showing pair energy plots ...\n")
  
  for tagk, tagj in ('min', 'min'), ('min', 'max'), ('max', 'max'):
    outPairEnergyGraph = f"{outPrefix}-epair-{tagk}-{tagj}.txt"
    bash(f"./plot_pair_energy.sh {outPairEnergyGraph}")
  
  if file_OK(outBegEnergyImage):
    bash(f"display -title '{outEndEnergyImage}' -filter Box -resize '{resize}%' {outEndEnergyImage}")

  Er(f"=== showing 2D slice of initial total energy {outBegTotalEnergyPlot} ...\n")
  bash(f"./plot_total_energy.sh {outBegTotalEnergyPlot}")

  Er(f"=== showing 2D slice of final total energy {outEndTotalEnergyPlot} ...\n")
  bash(f"./plot_total_energy.sh {outEndTotalEnergyPlot}")
  
  Er(f"=== making optimization movie ...\n")
  bash(f"./make_spots_movie.sh {outDir} {imageName}")
  
  exportExt = ("ppm" if NY == 400 else "pgm")
  exportName = f"{imageName}.{exportExt}"
  exportDir = f"{TEXTURE_BANK}/{exportExt}-{NX:03d}x{NY:03d}"
  exportFile = f"{exportDir}/{exportName}"

  Er(f"=== converting {outEndSpotsImage} --> {exportFile} ...\n")
  bash(f"convert {outEndSpotsImage} {exportFile}")
  bash(f"display -title 'projects: {exportName}' -filter Box -resize '{resize}%' {exportFile}")
  
  # ......................................................................

def single_clean(outPrefix):
  pat = f"{outPrefix}*." + "{png,txt,p?m}"
  Er(f"=== cleaning with 'rm -rf {pat}' ...\n")
  bash("rm -rf " + pat)
  # ......................................................................

def file_OK(path):
  return path != None and os.path.exists(path) and os.path.getsize(path) > 0
  # ......................................................................

def Er(msg):
  sys.stderr.write(msg)
  # ......................................................................

def run_command(command):
  result = subprocess.run(command, text = True, timeout = None)
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
