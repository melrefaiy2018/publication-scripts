c () { awk "BEGIN{ pi = 4.0*atan2(1.0,1.0); degree = pi/180.0; print $* }" ;}

if [ -z $1 ]
then
        echo "Please specify simulation time in ns. e.g. 6_md.sh 5" 
        exit
fi

ns=$1
steps=`c $ns*1000/0.002`
ps=`c $ns*1000`

module load GROMACS/2019-fosscuda-2018b 

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
neededFiles="${NEEDED_FILES:-$BASE_DIR}"
mdp_dir="${MDP_DIR:-$neededFiles/gromacs_mdp}"

cat "$mdp_dir/md.mdp" | sed s/XXXXXXX/$steps/g | sed s/YYYYY/$ps/ | sed s/ZZ/$ns/ > md.mdp

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_00_$ns.tpr

cat "$SCRIPT_DIR/gpujob.sh" | sed 's/XNAME/md/' > md.sh
echo "gmx mdrun -pin on -ntmpi 2  -ntomp 6 -nb gpu -pme cpu -bonded cpu -v -deffnm md_00_$ns" >> md.sh

sbatch md.sh
