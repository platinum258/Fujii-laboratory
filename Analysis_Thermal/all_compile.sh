
for no in 0 10
do
     sh clean.sh

     dev=$no sh setting.sh
     
     cd Program 
     sh clean.sh
     dev=$no sh compile.sh
     
     cd ..

      if [ -e ./Program/a.out ] ; then
           sh clean.sh
      else
           sh clean.sh
           echo '======================================================================='
           echo "dev= $no"
           echo '======================================================================='
           exit 1
      fi 

done

echo '======================================================================='
echo 'complete ' 
echo '======================================================================='

