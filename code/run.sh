#!/bin/bash

qcc -Wall -O2 program_v5.c -o program2 -L$BASILISK/gl -lglutils -lfb_tiny -lm 


rm data_output/*
#           ALPHA:  BETA:   Q:  M:  D_a:    c_0:    a_0:
./program2   1       2       1   1   1       2       1
#./program

python3 volume_grapher.py
#bash post_processing.sh

#python3 whole_system_animator.py
#ffmpeg -y -framerate 96 -i data_output/img%03d.png animation-single.mp4

mp4_file="animation-single.mp4"   
bash system_animator.sh $mp4_file 

rm program2


