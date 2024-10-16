#!/bin/bash
logfiles=$(ls data_output -1 | wc -l)
dt=$(ls data_output | sed 's/time=//' | sed 's/.dat//' | sed '2q;d' )
# use -persist to keep gnuplot up
gnuplot  <<-EOFMarker
    dt = ${dt}
    set terminal pngcairo size 512,512
    set size square
    set xrange [0:1]
    set yrange [0:2]
    
    do for [ii=0:${logfiles}] {
        
        infile = sprintf('data_output/time=%.3f.dat',ii*dt)
        outfile = sprintf('data_output/img%03d.png',ii)
        title_string = sprintf("time = %.3f",ii*dt)
        set output outfile
        
        set title title_string
        plot infile using 1:2 with filledcurve y1=0 title "c" , infile using 3:4 with filledcurve y1=0 title "a" 
        }
    
EOFMarker
# rest of script, after gnuplot exits

ffmpeg -y -framerate 48 -i data_output/img%03d.png $1
#rm data_output/*
