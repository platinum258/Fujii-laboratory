
#sh clean.sh
rm diff_tmp.sh list_diff.dat

#---
echo 'rm list_diff.dat 280.dat' > diff_tmp.sh
echo 'for fn in \' >> diff_tmp.sh

cd ./Program/
ls *.f90 >> ../diff_tmp.sh
cd ../

grep -l '.f90' ./diff_tmp.sh | xargs sed -i.bak -e 's/.f90/.f90 \\/g'

echo ' ' >> diff_tmp.sh
echo 'do' >> diff_tmp.sh
echo '   diff ${pwd1}/Program/${fn} $pwd2/Program/${fn} > 280.dat' >> diff_tmp.sh 
   
echo '   if [ -s 280.dat ]; then' >> diff_tmp.sh
echo '      echo "${fn}" >> list_diff.dat' >> diff_tmp.sh
echo '      rm 280.dat' >> diff_tmp.sh
echo '   else' >> diff_tmp.sh
echo '      rm 280.dat' >> diff_tmp.sh
echo '   fi' >> diff_tmp.sh
echo 'done' >> diff_tmp.sh

pwd1=${pwd1} pwd2=${pwd2} sh diff_tmp.sh 

rm diff_tmp.sh diff_tmp.sh.bak
#---
echo 'rm 280.dat' > diff_tmp.sh
echo 'for fn in \' >> diff_tmp.sh

cd ./Program/Parameters/
ls *.f90 >> ../../diff_tmp.sh
cd ../../

grep -l '.f90' ./diff_tmp.sh | xargs sed -i.bak -e 's/.f90/.f90 \\/g'

echo ' ' >> diff_tmp.sh
echo 'do' >> diff_tmp.sh
echo '   diff ${pwd1}/Program/Parameters/${fn} $pwd2/Program/Parameters/${fn} > 280.dat' >> diff_tmp.sh 
   
echo '   if [ -s 280.dat ]; then' >> diff_tmp.sh
echo '      echo "${fn}" >> list_diff.dat' >> diff_tmp.sh
echo '      rm 280.dat' >> diff_tmp.sh
echo '   else' >> diff_tmp.sh
echo '      rm 280.dat' >> diff_tmp.sh
echo '   fi' >> diff_tmp.sh
echo 'done' >> diff_tmp.sh

pwd1=${pwd1} pwd2=${pwd2} sh diff_tmp.sh 

rm diff_tmp.sh diff_tmp.sh.bak
#---


if [ ! -e ./list_diff.dat ] ;then
   echo '  '
   echo '================================================'
   echo 'No difference between all *.f90 in two project'
   echo '================================================'
   echo '  '
fi
