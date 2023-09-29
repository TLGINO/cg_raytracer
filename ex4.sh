#!/bin/bash

g++ main.cpp -o a.out

for i in $(seq 1 180); do
  ./a.out $i
done

input_pattern="%01d.ppm" # Change this pattern to match your file naming

output_file="output.mp4"
ffmpeg -framerate 60 -i "$input_pattern" -c:v vp9  -r 60 "$output_file" -y

rm [0-9]*.ppm
