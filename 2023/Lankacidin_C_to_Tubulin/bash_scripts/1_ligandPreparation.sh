
if [ -z $1 ]
then
	echo "Provide input file. ./1_ligandPreparation.sh complexFile"
	exit
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
neededFiles="${NEEDED_FILES:-$BASE_DIR}"
param_dir="${PARAM_DIR:-$neededFiles/ligand_paramatrization}"
ACPYPE_ROOT="${ACPYPE_ROOT:-$param_dir}"

module load  AmberTools/17-intel-2017b 

module load OpenBabel/2.4.0-intel-2016b-Python-2.7.12
cat $1 | grep "ATOM" | grep -v "LIG\|gdp\|GDP" > rec.pdb
cat $1 | grep "gdp\|GDP" > gdp.pdb
cat $1 | grep " LIG " > LIG.pdb


echo "###########################"
echo "###PARAMETERIZING LIGAND###"
echo "###########################"

mkdir LIG
cd LIG
cp ../LIG.pdb .


antechamber -i LIG.pdb -fi pdb -o LIG.mol2 -fo mol2 -c bcc -rn LIG -s 2
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod
rm A* leap.log sqm.* >& /dev/null

cp "$param_dir/LIGleap.in" .
tleap -f LIGleap.in

mkdir LIG.acpype
cp LIG_AC.* LIG.acpype
"$ACPYPE_ROOT/acpype.py" -i LIG.mol2 -b LIG
rm LIG.acpype/A* LIG.acpype/sqm.* >& /dev/null
cp LIG.acpype/LIG_GMX.* ..
cd ..

# processing gdp
echo ""
echo "##########################"
echo "####PARAMETERIZING GDP####"
echo "##########################"

mkdir gdp
cd gdp
cp ../gdp.pdb .
cp "$param_dir"/gdp* .
cp "$param_dir/frcmod.phos" .
cp "$param_dir/leaprc.ff99SBildn" .

tleap -f gdpleap.in

mkdir gdp.acpype
cp gdp_AC.* gdp.acpype
"$ACPYPE_ROOT/acpype.py" -i gdp.pdb -b gdp -n -3
rm gdp.acpype/A* gdp.acpype/sqm.* >& /dev/nul
cp gdp.acpype/gdp_GMX.* ..
cd ..
