#!/bin/bash

#SBATCH --job-name=matlab_thermalSZE_interacting_estates
#SBATCH --mem-per-cpu=15G
#SBATCH --partition=compute
#SBATCH --array=0-11%11
#SBATCH --mail-user=%u@oist.jp
#SBATCH --workdir=/work/BuschU/mickelli/lisa_code
#SBATCH --input=none
#SBATCH --output=/work/BuschU/mickelli/out/%j.out
#SBATCH --error=/work/BuschU/mickelli/error/%j.err
#SBATCH --time=4:00:00
## Load matlab (see 'module avail' for a list of available matlab versions)
module load matlab/R2016b

## matlab execution command
matlab_exe="matlab -nosplash -nodisplay -nojvm -nodesktop"

## parameterized matlab script execution command
#szarr=(1 2 3 4 5 6 7 8 9 10)
#g=(-1 0 1)
nb=(1 2 3 4 1 2 3 4 1 2 3 4)
g=(-1 -1 -1 -1 0 0 0 0 1 1 1 1)
gA=${g[$(($SLURM_ARRAY_TASK_ID))]}
#TA=${T[$(($SLURM_ARRAY_TASK_ID))]}
nbA=${nb[$(($SLURM_ARRAY_TASK_ID))]}

>&2 echo "[${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${gA}_${nbA} start]"
##ws_name="${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${msz}.mat"
mtalab_cmd="gAr=${gA};  nbAr=${nbA}; thermalSZE_InteractingBosons_estates_slurm"
#mtalab_cmd="n=${gA}; helloWorld"

## run matlabx
${matlab_exe} -r "${mtalab_cmd}"

## log completed job in .err file
>&2 echo "[${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${gA}_${nbA} completed]"
