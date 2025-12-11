module load GROMACS/2019-fosscuda-2018b 


neededFiles="/home/zc002u1/data/tubulin/neededFiles"

#restraining ligands

echo "0 & ! a H*
q" | gmx make_ndx -f LIG_GMX.gro -o index_LIG.ndx
echo "3" | gmx genrestr -f LIG_GMX.gro -n index_LIG.ndx -o posre_LIG.itp -fc 1000 1000 1000

echo "0 & ! a H*
q" | gmx make_ndx -f gdp_GMX.gro -o index_gdp.ndx
echo "3" | gmx genrestr -f gdp_GMX.gro -n index_gdp.ndx -o posre_gdp.itp -fc 1000 1000 1000

rest=`cat topol.top | grep "posre_LIG.itp" | wc -l`

if [ "$rest" -eq 0 ]
        then
                cat topol.top | sed 's/\[ moleculetype/; Ligand position restraints\n#ifdef POSRES\n#include "posre_LIG.itp"\n#endif\n\n\[ moleculetype/' > xx
                mv xx topol.top
		cat LIG_GMX.itp | sed 's/\[ moleculetype/; gdp position restraints\n#ifdef POSRES\n#include "posre_gdp.itp"\n#endif\n\n\[ moleculetype/' > yy
		mv yy LIG_GMX.itp


fi

echo "1 | 13 | 14 
q" | gmx make_ndx -f em.gro -o index.ndx

# NVT 
cp $neededFiles/nvt.mdp .
cp $neededFiles/npt.mdp .

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr

cat $neededFiles/gpujob.sh | sed 's/XNAME/equil/' > equil.sh
echo "gmx mdrun -pin on -ntomp 6 -v -deffnm nvt" >> equil.sh

echo "gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr" >> equil.sh

echo "gmx mdrun -pin on -ntomp 6 -v -deffnm npt" >> equil.sh

sbatch equil.sh


