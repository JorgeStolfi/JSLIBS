#! /usr/bin/python3 
# Last edited on 2025-04-09 05:18:34 by stolfi

from math import sin, cos, log, exp, pi, sqrt
import os, sys
import subprocess

PROG = "test_spots_image" 
TEXTURE_BANK = "projects/image-collections/texture-bank"

def Er(msg):
  sys.stderr.write(msg)
  # ......................................................................
  
def bash(cmd):
  # Execute the string {cmd} with "/bin/bash".
  subprocess.run([ cmd ], shell=True, executable="/bin/bash")
  # ......................................................................

def run_single(NX, NY, spotSize_min, spotSize_var, relSpotDist, maxSpots, numPasses):
  # Runs the program {PROG}.
  outDir = f"out/png-{NX:04d}x{NY:04d}"
  run_command([ "mkdir", "-pv", outDir ])

  outPrefix = f"{outDir}/spots-{spotSize_var}{spotSize_min}"
  imageName = f"spots-{spotSize_var}{spotSize_min}"
  
  outPairEnergyGraph = f"{outPrefix}-epair.txt"

  outBegTotalEnergyPlot = f"{outPrefix}-eplot-beg.txt"
  outEndTotalEnergyPlot = f"{outPrefix}-eplot-end.txt"

  outBegEnergyImage = f"{outPrefix}-espot-beg.png"
  outEndEnergyImage = f"{outPrefix}-espot-end.png"

  outBegSpotsImage = f"{outPrefix}-00000.png"
  outEndSpotsImage = f"{outPrefix}-final.png"

  Er(f"image size = {NX}x{NY} spot size = {spotSize_min} + {spotSize_var}\n")
  Er(f"relSpotDist = {relSpotDist:.3f}\n")
  
  single_clean(outPrefix)
  
  Er(f"=== making {outEndSpotsImage} ...\n")
  command = \
    [ PROG, 
      "-imageSize",    f"{NX}", f"{NY}",
      "-spotSize",     f"{spotSize_min}",  f"{spotSize_var}", 
      "-relSpotDist",  f"{relSpotDist:.4f}",
      "-maxSpots",     f"{maxSpots:d}",
      "-numPasses",    f"{numPasses:d}",
    ] 
  Er("command = [ " + " ".join(command) + " ]\n");
  run_command( command );

  run_command([ "ls", "-l", outBegSpotsImage, outEndSpotsImage, ]);
  
  resize = 200 if NY >= 200 else 400 # Percentage resize factor for display.

  Er(f"=== showing pair energy plot {outPairEnergyGraph} ...\n")
  run_command([ "./plot_pair_energy.sh", outPairEnergyGraph ])
  
  run_command([ "display", "-title", "'%f'", "-filter", "Box", "-resize", f"{resize}%", outBegEnergyImage, outEndEnergyImage, ])

  Er(f"=== showing 2D slice of initial total energy {outBegTotalEnergyPlot} ...\n")
  run_command([ "./plot_total_energy.sh", outBegTotalEnergyPlot ])

  Er(f"=== showing 2D slice of final total energy {outEndTotalEnergyPlot} ...\n")
  run_command([ "./plot_total_energy.sh", outEndTotalEnergyPlot ])
  
  Er(f"=== making optimization movie ...\n")
  run_command([ "./make_spots_movie.sh", outDir, imageName ])
  
  exportExt = ("ppm" if NY == 400 else "pgm")
  exportDir = f"{TEXTURE_BANK}/{exportExt}-{NX:03d}x{NY:03d}"
  exportFile = f"{exportDir}/{imageName}.{exportExt}"

  Er(f"=== converting {outEndSpotsImage} --> {exportFile} ...\n")
  run_command([ "convert", outEndSpotsImage, exportFile ])

  run_command([ "display", "-title", "'%f'", "-filter", "Box", "-resize", f"{resize}%", exportFile, outBegSpotsImage ])
  
  # ......................................................................

def run_command(command):
  result = subprocess.run(command, text = True)
  print(result.stderr)
  print(result.stdout)
  assert result.returncode == 0, f"** {command[0]} failed - returned status = {result}"
  return
  # ......................................................................

def single_clean(outPrefix):
  pat = f"{outPrefix}*." + "{png,txt,p?m}"
  Er(f"=== cleaning with 'rm -rf {pat}' ...\n")
  bash("rm -rf " + pat)
  # ......................................................................

def main():

  relSpotDist = 2.0

  # imageSize = 512
  # spotVars = range
  # numPasses = 100
  # for spotSize_min in spotSizes:
  #   for spotSize_var in range(0, 10-spotSize_min):
  #     run_single(imageSize, imageSize, spotSize_min, spotSize_var, relSpotDist, numPasses)
    
  # spotSize_min = 2
  
  spotSizes = ( (0,0), )
  imageSize = 512
  maxSpots = 9999
  numPasses = 100
  
  for spotSize_min,spotSize_var  in spotSizes:
    run_single(imageSize, imageSize, spotSize_min, spotSize_var, relSpotDist, maxSpots, numPasses)
  
  return
  # ......................................................................

main()
