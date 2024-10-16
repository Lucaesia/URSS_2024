#!/bin/bash

qcc -Wall -O2 program_v5.c -o program -L$BASILISK/gl -lglutils -lfb_tiny -lm 

#rm animation_output/*
#{start..increment..end}
for i in `seq 1.0 1.0 1.0`
do
    rm data_output/*
    #           ALPHA:  BETA:   Q:  M:  D_a:    c_0:    a_0
    ./program   1       1       1   1   1       1       1
    
    #python3 volume_grapher.py
    #bash post_processing.sh

    #python3 whole_system_animator.py
    #ffmpeg -y -framerate 24 -i data_output/img%03d.png animation.mp4
    mp4_file="animation_output/animation-M=${i}.mp4"
    
    bash system_animator.sh $mp4_file 
done
rm program


