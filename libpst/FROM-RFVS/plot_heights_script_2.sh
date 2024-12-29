heightmap="$1"; shift
texture="$1"; shift
png_output="$1"; shift

cp /home/saracchini/Raab/Versoes_2008_07_12/Fontes/plot_heights/plot_heights.pov plot_heights.pov
cp /home/saracchini/Raab/Versoes_2008_07_12/Fontes/plot_heights/camera.inc camera.inc
plot_heights ${heightmap} NONE ${texture} scene.inc
povray +A +R2 \
   Height=480 Width=640 \
   Library_Path=/home/saracchini/Raab/POVRAY/povlinux-3.6/povray-3.6/include \
   Input_File_Name=plot_heights.pov \
   Output_File_Name=${png_output}