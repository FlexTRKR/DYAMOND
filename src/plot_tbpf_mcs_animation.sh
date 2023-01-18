#!/bin/bash
# Makes Tb+PF MCS tracking quicklook animations over a subset domain and period

# Get config file from input
config=$1
# Get basename of config, remove suffix .yml
_config=$(basename "${config}" .yml)
# Subset strings to get the runname (e.g., config_dyamond_nicam.yml)
runname=${_config:15}
echo ${runname}
# Quicklook plot output directory
quicklook_dir=/global/cfs/cdirs/m1867/zfeng/dyamond-winter/quicklooks/${runname}/
# PyFLEXTRKR plotting code directory
# pyflex_dir=/global/homes/f/feng045/program/PyFLEXTRKR/Analysis/
# pyflex_dir='./'

# Domain boundary
figbasename='spcz_'
minlon=120.0
maxlon=230.0
minlat=-30.0
maxlat=15.0
# Start/end dates
startdate='2020-02-01T00'
enddate='2020-02-10T00'

source activate flextrkr

python plot_subset_tbpf_mcs_tracks_demo.py -s ${startdate} -e ${enddate} \
    -c ${config} \
    -o vertical --figsize 17 4 \
    -p 1 --nprocesses 32 \
    --extent ${minlon} ${maxlon} ${minlat} ${maxlat} \
    --output ${quicklook_dir} \
    --figbasename ${figbasename}

# echo 'Making animations from quicklook plots ...'
# ffmpeg -framerate 2 -pattern_type glob -i ${quicklook_dir}'*.png' \
#     -c:v libx264 -r 10 -crf 20 -pix_fmt yuv420p \
#     -y ${quicklook_dir}quicklook_animation_${runname}.mp4
# echo 'View animation here: '${quicklook_dir}quicklook_animation_${runname}.mp4