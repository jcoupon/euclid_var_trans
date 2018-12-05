#! /bin/bash
set -o nounset

# ----------------------------------------------------------------------------------- #
# options
# ----------------------------------------------------------------------------------- #

TMPDIR=$HOME/data/tmp
STILTS='java -Xmx2048M -jar /Users/coupon/local/bin/stilts.jar'

COSMOSDIR=$HOME/data/COSMOS/COSMOS2015
DATADIR=$HOME/data/euclid/varTrans
HOST=$( hostname -s )

# default options
MPIRUN=mpirun
SWOT=$HOME/local/source/GitHub/swot/bin/swot
VENICE=$HOME/local/source/GitHub/venice/bin/venice
NP=2

if [[ "$HOST" =~ "ojingo" ]]; then
   MPIRUN=/opt/openmpi-1.8.6_clang/bin/mpirun
   NP=8
fi
if [[ "$HOST" =~ "hayabusa" ]]; then
   MPIRUN=/opt/openmpi-1.8.5/bin/mpirun
   NP=4
fi

# WMAP9 {Omega_m, Omega_b, Omega_Lambda, sigma_8, n_s, h} = {0.2793, 0.0463, 0.7207, 0.821, 0.972, 0.700}
H0=70.0; OM=0.2793; OL=0.7207 # WMAP Hinshaw et al. (2009)
