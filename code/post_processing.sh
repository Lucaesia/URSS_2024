#!/bin/bash
minval=0    # the result of some (omitted) calculation
maxval=4219   # ditto

# use -persist to keep gnuplot up
gnuplot  <<-EOFMarker
    dt = 0.002
    set terminal pngcairo size 512,512
    set size square
    set xrange [0:1]
    set yrange [-1:1.1]
    
    do for [ii=0:500] {
        
        infile = sprintf('data_output/time=%.3f.dat',ii*dt)
        outfile = sprintf('data_output/img%03d.png',ii)
        title_string = sprintf("time = %.3f",ii*dt)
        set output outfile
        
        set title title_string
        plot infile using 1:3 with lines title "a" 
        }
    
EOFMarker
# rest of script, after gnuplot exits

ffmpeg -y -framerate 96 -i data_output/img%03d.png output_a.mp4 
#rm data_output/*
