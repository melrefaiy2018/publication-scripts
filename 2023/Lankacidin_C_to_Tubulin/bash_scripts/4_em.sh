module load GROMACS/2019-foss-2018b 

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
neededFiles="${NEEDED_FILES:-$BASE_DIR}"
mdp_dir="${MDP_DIR:-$neededFiles/gromacs_mdp}"

cp "$mdp_dir/em.mdp" .
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

cat "$SCRIPT_DIR/gpujob.sh" | sed 's/XNAME/em/' > em.sh
echo "gmx mdrun -pin on -ntomp 1 -v -deffnm em" >> em.sh

sbatch em.sh
