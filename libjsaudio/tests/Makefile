# Last edited on 2023-03-02 13:38:47 by stolfi

IGNOREDIRS :=
  
include ${STOLFIHOME}/programs/c/GENERIC-ROOT-DIR.make
  
######################################################################
# OGG to AU conversion:

in/%-raw.au: in/%-raw.ogg 
	sox -V in/$*-raw.ogg -s -w -r 44100 in/$*-raw.au polyphase

