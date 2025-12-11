c () { awk "BEGIN{ pi = 4.0*atan2(1.0,1.0); degree = pi/180.0; print $* }" ;}

if [ -z $1 ] || [ -z $2 ]
then
        echo "Please specify input file and simulation time in ns. e.g. 6_md.sh previousXtc 5" 
        exit
fi

ns=$2
steps=`c $ns*1000/0.002`
ps=`c $ns*1000`

name=`basename $1 .xtc`
st=`echo $1 | cut -d"_" -f3 | cut -d"." -f1`
et=`c "$st+$ns"`

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
neededFiles="${NEEDED_FILES:-$BASE_DIR}"
mdp_dir="${MDP_DIR:-$neededFiles/gromacs_mdp}"

cat "$mdp_dir/md.mdp" | sed s/XXXXXXX/$steps/g | sed s/YYYYY/$ps/ | sed s/ZZ/$ns/ > md.mdp

gmx grompp -f md.mdp -c $name.gro -t $name.cpt -p topol.top -n index.ndx -o md_"$st"_"$et".tpr

cat "$SCRIPT_DIR/gpujob.sh" | sed 's/XNAME/md/' > md.sh
echo "gmx mdrun -pin on -ntmpi 2 -ntomp 6 -nb gpu -pme cpu -bonded cpu -v -deffnm md_"$st"_"$et"" >> md.sh

sbatch md.sh
