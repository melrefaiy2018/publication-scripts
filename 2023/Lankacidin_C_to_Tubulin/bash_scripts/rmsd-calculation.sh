
gmx rms -s md_00_50.tpr -f md_00_50_center.xtc -o rmsd.xvg -tu ns
gmx gyrate -s md_00_50.tpr -f md_00_50_center.xtc -o gyrate.xvg
gmx rmsf -f md_00_50_center.xtc -s md_00_50.gro -n index.ndx -res -fit
