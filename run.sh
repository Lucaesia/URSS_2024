#!/bin/bash

qcc -Wall -O2 program_v1.c -o program -L$BASILISK/gl -lglutils -lfb_tiny -lm 


./program  > "boundary_pos.dat"
rm program
#bash post_processing.sh
