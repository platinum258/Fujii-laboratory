#!/bin/sh
echo "Hello, gifWorld!"
#echo "======================================="
#echo "              Start trim               "
#echo "======================================="
#for gnumber in *_Configuration.ps; do
#done
#cd Result_"${TODAY}" 
#for fname in *.ps; do
#  newfname=$(echo "${fname}" | sed "s/\.ps$/_trim.ps/g")
#  convert -trim "${fname}" "${newfname}"
#done
echo "======================================="
echo "         Start change to gif           "
echo "======================================="
for fname in *.ps; do
  newfname=$(echo "${fname}" | sed "s/\.ps$/.gif/g")
  convert "${fname}" "${newfname}"
done
echo "======================================="
echo "         Make gif animation            "
echo "======================================="
    mkdir Animation
#    for i in 'seq 15'
#     cp "${gnumber}"_"${i}"_Configuration_trim.gif
#     cp "${gnumber}"_"${i}"_1_1_Mises_Stress_trim.gif  
#    done 
    convert -adjoin *_Structure.gif Animation/Configuration.gif
    #convert -adjoin *_Mises_Stress_trim.gif  Animation/Distoribution.gif
    mkdir Image
    mv *.ps Image/
    mv *.gif Image/
echo "======================================="
echo "         END Make Animation!!          "
echo "======================================="
