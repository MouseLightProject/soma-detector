#!/bin/bash
cd /groups/mousebrainmicro/mousebrainmicro/scripts/soma-detector/
out=soma-detector.out.round1.txt
bsub -P mouselight -n 8 -oo $out -eo $out /misc/local/matlab-2019a/bin/matlab -nodisplay -batch "modpath; find_somata_2020_12_31_new_v1"

