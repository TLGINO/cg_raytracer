#!/bin/bash

g++ main.cpp -Ofast -o a.out




for i in $(seq 1 100); do
    result=$(echo "2 * s(sqrt($i))" | bc -l)
    formatted_result=$(printf "%.4f" $result)

    echo $formatted_result

    ./a.out $i $formatted_result $formatted_result
done




input_pattern="%01d.ppm" # Change this pattern to match your file naming

output_file="output.mp4"
ffmpeg -framerate 60 -i "$input_pattern" -c:v vp9  -r 60 "$output_file" -y

rm [0-9]*.ppm