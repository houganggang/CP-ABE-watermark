#!/bin/bash
module load matlab/2017a

mkdir -p compiledmatlab

cat > build.m <<END
addpath(genpath('/N/u/brlife/git/jsonlab'))  
addpath(genpath('/N/u/brlife/git/geom3d'))
mcc -m -R -nodisplay -d compiledmatlab main
exit
END
matlab -nodisplay -nosplash -r build
