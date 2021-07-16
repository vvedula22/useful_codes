#!/bin/bash
nstart=0
nend=6
nfreq=1
n1=20
nfiles=$(( $(($nend-$nstart)) / $nfreq + 1))
fileTemplate=$(printf "_registered")
for (( i=0; i < $nfiles; i++ ))
do
    ntime1=$(($nstart + $i * $nfreq))
    ntime2=$(($n1 + $i * $nfreq))
    fname1=$(printf "%02d" $ntime1)$fileTemplate$(printf ".vtk")
    fname2=$(printf "%02d" $ntime2)$fileTemplate$(printf ".vtk")
    echo $fname1     $fname2
    mv $fname1  $fname2
done


