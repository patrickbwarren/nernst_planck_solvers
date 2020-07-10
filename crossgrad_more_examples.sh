# Fluorescein (disodium salt) gradient crossed with KCl gradient

./crossgrad_creator.py --ionx=Na2Fl1 --iony=K1Cl1 --nside=20 --nsave=5 --nstep=250 fluorsc --run
./crossgrad_solver fluorsc.rc 
ls fluorsc_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=fluorsc.rc --save 
ffmpeg  -i fluorsc_f%05d.dat.png -c:v libx264 -r 30 -pix_fmt yuv420p fluorsc.mp4

# NaCl gradient crossed with KCl gradient

./crossgrad_creator.py --ionx=Na1Cl1 --iony=K1Cl1 --nside=20 --nsave=5 --nstep=250 nakcl --run
./crossgrad_solver nakcl.rc 
ls nakcl_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=nakcl.rc --save 
ffmpeg  -i nakcl_f%05d.dat.png -c:v libx264 -r 30 -pix_fmt yuv420p nakcl.mp4

# NaOH gradient crossed with KBr gradient

./crossgrad_creator.py --ionx=Na1OH1 --iony=K1Br1 --nside=20 --nsave=5 --nstep=250 --dt=0.025 naohkbr --run
./crossgrad_solver naohkbr.rc 
ls naohkbr_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=naohkbr.rc --save 
ffmpeg  -i naohkbr_f%05d.dat.png -c:v libx264 -r 30 -pix_fmt yuv420p naohkbr.mp4

#./crossgrad_creator.py --ionx="{'Na':1, 'Cl':1}" --iony="{'K':1, 'Cl':1}" --ionb="{'Na':2, 'Fl':1}" --back=0.01 --nside=20 --nsave=10 tracer
#./crossgrad_solver tracer.rc 
#ls tracer_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=tracer.rc --save

# LiCl gradient crossed with PbCO3 gradient

./crossgrad_creator.py --ionx=Li1Cl1 --iony=Pb1COt1 --nside=20 --nsave=5 --nstep=250 liclpbco3 --run
./crossgrad_solver liclpbco3.rc 
ls liclpbco3_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=liclpbco3.rc --save 

# MAIN EXAMPLES USED IN PAPER  <====  USE THIS EXAMPLE !

# HCl gradient crossed with KBr gradient -- MEDIUM (20)

./crossgrad_creator.py --ionx=H1Cl1 --iony=K1Br1 --nside=20 --nsave=8 --nstep=800 --dt=0.025 hclkbrM --run
./crossgrad_solver hclkbrM.rc 
ls hclkbrM_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=hclkbrM.rc --save 
ffmpeg  -i hclkbrM_f%05d.dat.png -c:v libx264 -r 30 -pix_fmt yuv420p hclkbrM.mp4

# HCl gradient crossed with KBr gradient -- LARGE (40) 

./crossgrad_creator.py --ionx=H1Cl1 --iony=K1Br1 --nside=40 --nsave=32 --nstep=3200 --dt=0.025 hclkbrL --run
./crossgrad_solver hclkbrL.rc 
ls hclkbrL_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=hclkbrL.rc --save 
ffmpeg  -i hclkbrL_f%05d.dat.png -c:v libx264 -r 30 -pix_fmt yuv420p hclkbrL.mp4

# HCl gradient crossed with KBr gradient -- EXTRA LARGE (50) - 56 minute run time

./crossgrad_creator.py --ionx=H1Cl1 --iony=K1Br1 --nside=50 --nsave=50 --nstep=5000 --dt=0.025 hclkbrXL --run
./crossgrad_solver hclkbrXL.rc 
ls hclkbrXL_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=hclkbrXL.rc --save 
ffmpeg  -i hclkbrXL_f%05d.dat.png -c:v libx264 -r 30 -pix_fmt yuv420p hclkbrXL.mp4

# HCl gradient crossed with KBr gradient -- EXTRA EXTRA LARGE (100) - 38 hour run time

./crossgrad_creator.py --ionx=H1Cl1 --iony=K1Br1 --nside=100 --nsave=125 --nstep=800*25 --dt=0.025 hclkbrXXL --run
./crossgrad_solver hclkbrXXL.rc 
ls hclkbrXXL_f*.dat | xargs -n 1 ./crossgrad_viewer.py --rc=hclkbrXXL.rc --save 
ffmpeg  -i hclkbrXXL_f%05d.dat.png -c:v libx264 -r 30 -pix_fmt yuv420p hclkbrXXL.mp4

# post processing to extract phoretic trajectories

./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-60 | gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' \
    | sort -g -k3 -k2 > hclkbrXXL_n60_all_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-40 | gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' \
    | sort -g -k3 -k2 > hclkbrXXL_n40_all_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-20 | gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' \
    | sort -g -k3 -k2 > hclkbrXXL_n20_all_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=20 | gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' \
    | sort -g -k3 -k2 > hclkbrXXL_p20_all_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=40 | gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' \
    | sort -g -k3 -k2 > hclkbrXXL_p40_all_traj.dat

./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-60 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_n60_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-40 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_n40_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-20 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_n20_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=20 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_p20_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=40 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_p40_50-50_single_traj.dat

./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-60 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_n60_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-40 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_n40_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=-20 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_n20_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=20 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_p20_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrXXL.rc --zeta=40 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXXL_p40_20-20_single_traj.dat

./crossgrad_trajectory.py hclkbrXL.rc --zeta=-60 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_n60_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrXL.rc --zeta=-40 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_n40_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrXL.rc --zeta=-20 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_n20_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrXL.rc --zeta=20 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_p20_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrXL.rc --zeta=40 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_p40_50-50_single_traj.dat

./crossgrad_trajectory.py hclkbrXL.rc --zeta=-60 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_n60_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrXL.rc --zeta=-40 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_n40_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrXL.rc --zeta=-20 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_n20_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrXL.rc --zeta=20 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_p20_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrXL.rc --zeta=40 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrXL_p40_20-20_single_traj.dat

./crossgrad_trajectory.py hclkbrM.rc --zeta=-60 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_n60_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrM.rc --zeta=-40 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_n40_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrM.rc --zeta=-20 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_n20_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrM.rc --zeta=20 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_p20_50-50_single_traj.dat
./crossgrad_trajectory.py hclkbrM.rc --zeta=40 --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_p40_50-50_single_traj.dat

./crossgrad_trajectory.py hclkbrM.rc --zeta=-60 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_n60_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrM.rc --zeta=-40 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_n40_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrM.rc --zeta=-20 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_n20_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrM.rc --zeta=20 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_p20_20-20_single_traj.dat
./crossgrad_trajectory.py hclkbrM.rc --zeta=40 --xy='[20,20]' --single | \
    gawk '/position/ { printf "%g\t%g\t%g\t%g\t%g\n", $1, $2, $3, $8, $9 }' > hclkbrM_p40_20-20_single_traj.dat


