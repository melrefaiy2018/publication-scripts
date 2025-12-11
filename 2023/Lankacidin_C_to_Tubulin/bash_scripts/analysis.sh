module load GROMACS/2019-foss-2018b 



if [ -z $1 ] && [ -z $2 ] && [ -z $3 ]
then
	echo "Provide the trajectroy file, followed by the .gro file, then em.tpr"
	exit
fi
mkdir analysis/

echo "#######################"
echo "####Receptor###########"
echo "#######################"

echo "3
3" | gmx rms -s $3 -f $1 -o rmsd-rec.xvg -tu ns

echo "1"| gmx gyrate -s $3 -f $1 -o gyrate-rec.xvg
echo "3"|gmx rmsf -f $1 -s $2 -n index.ndx -res -fit -o rmsf-rec.xvg

echo "#######################"
echo "########Ligand#########"
echo "#######################"

echo "13
13" | gmx rms -s $3 -f $1 -o rmsd-LIG.xvg -tu ns
echo "13" |gmx gyrate -s $3 -f $1 -o gyrate-LIG.xvg


mv rmsd-rec.xvg gyrate-rec.xvg rmsd-LIG.xvg  gyrate-LIG.xvg rmsf-rec.xvg  analysis/




