#!/usr/bin/bash

# NaCl disc overlapped with a KBr disc - this generates sample commands to do a run

jobname=xdisclarge

size=200
shape=005,0.125,0.025
run_control="--njacobi=10000 --niter=40"
grad=0.1,1

echo ./crossgrad_simple.py $jobname --solve --nx=$size --ny=$size --shape=$shape --cx=$grad --cy=$grad $run_control

echo ./crossgrad_viewer.py --rc ${jobname}.rc --show --phi-max=10.0 --detail --zoom 1.2
