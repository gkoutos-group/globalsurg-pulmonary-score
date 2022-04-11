#!/bin/bash
#SBATCH --account=gkoutosg-variant-prediction
#SBATCH --qos castles
#SBATCH --partition=icelake-castles
#SBATCH --mail-type ALL
#SBATCH --nodes 1
#SBATCH --ntasks 70
#SBATCH --time 6:0:0
#SBATCH --mem 100G
#SBATCH --job-name Ice

module purge;
module load bluebear
module load R/4.0.0-foss-2020a

Rscript /rds/projects/g/gkoutosg-variant-prediction/VictorSaisa/FinalCode/March11_PulmonaryComplications.R
