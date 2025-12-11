

et=`ls md*xtc | tail -1 | cut -d"_" -f3 | cut -d"." -f1`
gmx trjcat -f md_*xtc -o md_00_"$et".xtc --settime
