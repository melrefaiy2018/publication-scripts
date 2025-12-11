module load GROMACS/2019-foss-2018b 


neededFiles="/home/zc002u1/data/tubulin/neededFiles"

cp $neededFiles/em.mdp .
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

cat $neededFiles/gpujob.sh | sed 's/XNAME/em/' > em.sh
echo "gmx mdrun -pin on -ntomp 1 -v -deffnm em" >> em.sh

sbatch em.sh

