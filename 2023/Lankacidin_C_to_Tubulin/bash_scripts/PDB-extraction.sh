if [ -z $1 ] && [ -z $2 ]
then
	echo "Provide trajectory file, follwed by time(10000)"
	exit
fi
module load GROMACS/2019-fosscuda-2018b 
mkdir PDB/
gmx trjconv -s em.tpr -f $1 -n index2.ndx -sep -dt $2 -o PDB/PDB.pdb 
