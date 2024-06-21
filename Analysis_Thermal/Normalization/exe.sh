

dev=$dev parameter_file=$parameter_file sh compile.sh

if [ ! -e a.out ] ; then
      echo '======================================================================='
      echo "No a.out : Directory $PWD"
      echo 'exe.sh 8'
      echo '======================================================================='

      exit 1
      exit 
fi

#if [ -e ./a.out ] ; then
if [ -e a.out ] ; then
     ./a.out
else
      echo '======================================================================='
      echo "No a.out : Directory $PWD"
      echo 'exe.sh 19'
      echo '======================================================================='

      exit 1
fi 

#gprof ./a.out gmon.out > gmon.log
rm a.out

#sh diff.sh

if [ -e ./Result*.dat ] ; then
      echo "Calculation FINISHED : Directory : $PWD"
else
      echo '======================================================================='
      echo "No Result.dat : Directory : $PWD"
      echo '======================================================================='

      exit 1
fi

