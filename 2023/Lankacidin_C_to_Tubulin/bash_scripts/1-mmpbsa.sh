if [ -z $1 ]
then
	echo "Provide input xtc file."
	exit
fi

name=`basename $1 .xtc`
module load GROMACS/2016.3-foss-2016b-GPU-enabled

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
neededFiles="${NEEDED_FILES:-$BASE_DIR}"
mdp_dir="${MDP_DIR:-$neededFiles/gromacs_mdp}"

echo "1 | 14
q" | gmx make_ndx -f em.gro -n index.ndx -o index2.ndx

mkdir mmpbsa
cd mmpbsa
cp "$SCRIPT_DIR/2-mmpbsa.sh" .


echo "1
0" | gmx trjconv -s ../em.tpr -f ../$1 -o ../"$name"_center.xtc -center -pbc mol -ur compact

cp "$mdp_dir/pbsa.mdp" .
cp "$SCRIPT_DIR/mmpbsajob.sh" .
