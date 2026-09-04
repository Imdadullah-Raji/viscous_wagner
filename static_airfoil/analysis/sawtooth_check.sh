#!/bin/bash
# Detect the 2*dt sawtooth described in FINDINGS.md Part 4.
# ~100% step-to-step reversals in Cl means the run is corrupted even though it
# completed normally. Usage: sawtooth_check.sh <case-dir> [t_min]
#
# Concatenates ALL forceCoeffs restart dirs and sorts numerically by time.
# (Do not use `find | sort | tail -1`: the path ".../0.5/" sorts BEFORE ".../0/"
# because '.' < '/', so that silently picks the wrong file on a restarted run.)
d="$1"; w="${2:-10}"
files=$(find "$d/postProcessing/forceCoeffs" -name '*.dat' 2>/dev/null)
[ -z "$files" ] && { echo "  no forceCoeffs.dat under $d"; exit 1; }
cat $files | grep -v '^#' | sort -g -k1 | awk '!seen[$1]++' | awk -v w="$w" '
$1>w {cl[n++]=$4}
END{
  if(n<50){print "  insufficient data (n="n")"; exit}
  f=0; for(i=2;i<n;i++){d1=cl[i]-cl[i-1]; d2=cl[i-1]-cl[i-2]; if(d1*d2<0) f++}
  mn=cl[0];mx=cl[0];s=0
  for(i=0;i<n;i++){s+=cl[i]; if(cl[i]<mn)mn=cl[i]; if(cl[i]>mx)mx=cl[i]}
  r=100*f/(n-2)
  printf "  n=%d  Cl_mean=%.4f  range=[%.4f, %.4f]  p-p=%.2e  step-reversals=%.1f%%  %s\n",
         n, s/n, mn, mx, mx-mn, r, (r>80 ? "<-- SAWTOOTH, RESULT IS GARBAGE" : "ok")
}'
