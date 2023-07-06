#!/usr/bin/bash

# NaCl disc overlapped with a KBr disc or Gaussian overlapped with a Gaussian
# - this generates sample commands to do a run.

size=200
shape=0.05,0.125,0.025
run_control="--njacobi=10000 --niter=40"
grad=0.1,1

jobname=xdisclarge
echo ./crossgrad_simple.py $jobname --solve --nx=$size --ny=$size --shape=$shape --cx=$grad --cy=$grad $run_control
echo ./crossgrad_viewer.py --rc ${jobname}.rc --show --phi-max=10.0 --detail --zoom 1.2

shape=0.05,0.25,0.025
jobname=xgausslarge
echo ./crossgrad_simple.py $jobname --solve --nx=$size --ny=$size --shape=$shape --cx=$grad --cy=$grad $run_control --gauss
echo ./crossgrad_viewer.py --rc ${jobname}.rc --show --phi-max=10.0 --detail --zoom 1.2
