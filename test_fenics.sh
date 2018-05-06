#!/bin/bash
# Job name:
#SBATCH --job-name=rmcl5_9
##ro2l2s9
###rp25o2s9+
#
# Project:
#SBATCH --account=nn9316k
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
# Max memory usage per task:
#SBATCH --mem-per-cpu=24G
#
# Number of tasks (cores):
#SBATCH --ntasks=1



## Setup job environment:
source /cluster/bin/jobsetup
module purge    # clear any inherited modules
set -o errexit  # exit on errors
source /usit/abel/u1/johannr/fenics-2016.2.0.abel.gnu.conf
##source /usit/abel/u1/johannr/fenics-2016.1.0.abel.gnu.conf
##source /usit/abel/u1/johannr/fenics-1.6.0.abel.gnu.conf
##source /usit/abel/u1/johannr/fenics-1.5.0.abel.gnu.conf

## Set up input and output files:
cp $SUBMITDIR/*.xml $SCRATCH
cp $SUBMITDIR/*.py $SCRATCH

chkfile "UQ_*"
chkfile "displacement*"
chkfile "subdomains*"


cd $SCRATCH
#python demo-dolfin.py
#python main_simple.py
#python model_onlyparams.py
#python  uqanalysis_pce_rfield.py
python uqanalysis_mc_rfield_bayer_mod.py
#python complexmodel.py
#python complexmodel_pce-pcm.py
