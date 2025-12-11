#! /bin/bash
#SBATCH --job-name=mmpbsa
#SBATCH --partition=cpu
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=s-mohamed.elrefaiy@zewailcity.edu.eg
#SBATCH --mail-type=END
#SBATCH --time=16:00:00

module load GROMACS/2019-fosscuda-2018b

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
bfe_dir="${BFE_DIR:-$BASE_DIR/binding_free_energy}"
PBSA_BIN="${PBSA_BIN:-$bfe_dir/g_mmpbsa/bin/g_mmpbsa}"


gmx grompp -f ../md.mdp -c ../em.gro -n ../index.ndx -p ../topol.top -o topol.tpr

cp ../index2.ndx .


echo "27
13" | $PBSA_BIN -f ../md_00_50_center.xtc -s topol.tpr -n index2.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp -dt 250 -b 5000
