#! /bin/bash
#SBATCH --job-name=mmpbsa
#SBATCH --partition=cpu
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=s-mohamed.elrefaiy@zewailcity.edu.eg
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --time=16:00:00

module load Python/2.7.12-foss-2016b 
module load GROMACS/2016.3-foss-2016b-GPU-enabled
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
bfe_dir="${BFE_DIR:-$BASE_DIR/binding_free_energy}"
PBSA_BIN="${PBSA_BIN:-$bfe_dir/g_mmpbsa/bin/g_mmpbsa}"



echo "26
14 | $PBSA_BIN -f ../md_00_50_center.xtc -s ../em.tpr -n ../index2.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp -dt 250 " >> mmpbsajob.sh

sbatch mmpbsajob.sh

~                                                        
