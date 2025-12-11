#! /bin/bash
#SBATCH --job-name=XNAME
#SBATCH --partition=gpu
#SBATCH --account=g.zc002
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mail-user=s-mohamed.elrefaiy@zewailcity.edu.eg
#SBATCH --mail-type=END
#SBATCH --time=16:00:00

module load GROMACS/2019-fosscuda-2018b

