#!/bin/bash

min=$1
max=$2
in_1=$3
in_2=$4
out=$5

paste $in_1 $in_2 | \
awk '{print $1, $2, $5, $11}' | \
awk -F '[:]' '{print $1, $2, $3, $4}' | \
awk -v min="$min" -v max="$max" '{ if ($4 >= max && $6 <= min || $4 <= min && $6 >= max) print $1"\t"$2}' | \
awk '/LG/' > $out