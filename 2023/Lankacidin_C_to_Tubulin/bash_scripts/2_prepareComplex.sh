# prepare complex and topology

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
neededFiles="${NEEDED_FILES:-$BASE_DIR}"
param_dir="${PARAM_DIR:-$neededFiles/ligand_paramatrization}"

cp gdp/gdp.acpype/gdp_GMX.gro .
cp "$param_dir/gdp_GMX.itp" .
cp LIG/LIG.acpype/LIG_GMX.gro .
cp LIG/LIG.acpype/LIG_GMX.itp .

module load GROMACS/2019-fosscuda-2018b 



rm \#* topol.top rec_processed.gro >& /dev/null

echo "6
1" | gmx pdb2gmx -f rec.pdb -o rec_processed.gro -ignh

cat LIG_GMX.gro gdp_GMX.gro | grep " 1 " > LIGgdp.tmp
lignr=`cat LIGgdp.tmp | wc -l`
echo "Good ROcking Metal Altar for Chronical Sinners" > complex.gro
cat rec_processed.gro | head -2 | tail -1 | awk -v n=$lignr '{print " "$1+n}' >> complex.gro

head -n -1 rec_processed.gro | tail -n +3 >> complex.gro
cat LIGgdp.tmp >> complex.gro
cat rec_processed.gro | tail -1 >> complex.gro

rm LIGgdp.tmp

cat topol.top | sed 's/forcefield.itp\"/forcefield.itp\"\n\n; Include ligand topology\n#include "LIG_GMX.itp"/' > x
mv x topol.top
echo "LIG                 1" >> topol.top
echo "gdp                 1" >> topol.top

cat LIG_GMX.itp | sed 's/\[ moleculetype/; Include gdp topology\n#include "gdp_GMX.itp"\n\n\[ moleculetype/' > xx
mv xx LIG_GMX.itp
