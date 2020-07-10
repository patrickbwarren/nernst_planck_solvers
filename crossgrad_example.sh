#!/usr/bin/bash

# HCl gradient crossed with KBr gradient -- MEDIUM (20)

header=hclkbrMs

./crossgrad_creator.py --ionx=H1Cl1 --iony=K1Br1 --nside=20 --nsave=8 --nstep=800 --steady --dt=0.025 ${header} --run
./crossgrad_solver ${header}.rc
ls hclkbrM_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=${header}.rc --save
ffmpeg  -i hclkbrM_f%05d.dat.png -c:v libx264 -r 30 -pix_fmt yuv420p ${header}.mp4
