
if [ -z $1 ]
then
	echo "Provide input xtc file."
	exit
fi

name=`basename $1 .xtc`
module load GROMACS/2019-fosscuda-2018b 


neededFiles="/home/zc002u1/data/tubulin/neededFiles"
pbsa="/home/zc002u1/data/tubulin/g_mmpbsa/bin"

echo "1 | 14
q" | gmx make_ndx -f em.gro -n index.ndx -o index2.ndx

mkdir mmpbsa
cd mmpbsa

echo "1
0" | gmx trjconv -s ../em.tpr -f ../$1 -o ../"$name"_center.xtc -center -pbc mol -ur compact

cp $neededFiles/pbsa.mdp .
cp $neededFiles/mmpbsajob.sh .

