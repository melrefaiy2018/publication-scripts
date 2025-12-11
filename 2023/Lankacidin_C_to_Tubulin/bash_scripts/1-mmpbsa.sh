if [ -z $1 ]
then
	echo "Provide input xtc file."
	exit
fi

name=`basename $1 .xtc`
module load GROMACS/2016.3-foss-2016b-GPU-enabled

neededFiles="/home/zewail002u1/data/tubulin/neededFiles"
pbsa="/home/zewail002u1/data/tubulin/g_mmpbsa/bin"

echo "1 | 14
q" | gmx make_ndx -f em.gro -n index.ndx -o index2.ndx

mkdir mmpbsa
cd mmpbsa
cp $neededFiles/2-mmpbsa.sh .


echo "1
0" | gmx trjconv -s ../em.tpr -f ../$1 -o ../"$name"_center.xtc -center -pbc mol -ur compact

cp $neededFiles/pbsa.mdp .
cp $neededFiles/mmpbsajob.sh .

