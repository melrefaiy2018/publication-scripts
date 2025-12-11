module load GROMACS/2019-fosscuda-2018b 


neededFiles="/home/zc002u1/data/tubulin/neededFiles"

gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

echo "#################TOP#################"
tail topol.top
echo "#####################################"

cp $neededFiles/ions.mdp .
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo 16 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

echo "#################TOP#################"
tail topol.top
echo "#####################################"

