#!/bin/bash

# extract_lift_coeffs.sh [caseRoot] [out.csv]  (default: NACA0012/steady_re500 -> lift_coeffs.csv)
src="${1:-NACA0012/steady_re500}"
out="${2:-lift_coeffs.csv}"
> "$out"

for dir in "$src"/steady_aoa*/; do
    file="${dir}postProcessing/forceCoeffs/0/forceCoeffs.dat"
    
    [[ -f "$file" ]] || continue

    aoa="${dir%/}"
    aoa="${aoa##*/steady_aoa}"

    cl=$(tail -n 1 "$file" | awk '{print $4}')

    printf "%s %s\n" "$aoa" "$cl" >> "$out"
done

sort -g -o "$out" "$out"
