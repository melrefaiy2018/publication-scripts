cat polar.xvg | grep -v "@\|#" | awk '{E=$4-$3-$2} {print $1"   "E}' > pb
cat energy_MM.xvg | grep -v "@\|#" | awk '{print $1"   "$8}' > mm
cat energy_MM.xvg | grep -v "@\|#" | awk '{print $1"   "$7}' > ele
cat energy_MM.xvg | grep -v "@\|#" | awk '{print $1"   "$6}' > vdw
cat apolar.xvg | grep -v "@\|#" | awk '{E=$4-$3-$2}{print $1"    "E}' > sa
paste mm pb sa | awk '{E=$2+$4+$6}{print $1"    "E}' > sum
