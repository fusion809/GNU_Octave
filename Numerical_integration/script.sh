#!/bin/bash
start_time="$(date +%s)"
octave --eval 'run("Chebyshev-Gauss_quadrature_Gaussian.m")'
time_taken="$(($(date +%s)-start_time))"
echo "It took ${time_taken} seconds for the script to run"
