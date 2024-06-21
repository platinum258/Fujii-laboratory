
filename='Device_Number'

#===========================================================================
#if [ ./Parameters.f90 -ot ./Parameters/$parameter_file ] ; then
     dev=$dev parameter_file=$parameter_file sh update_Parameters.sh
#fi

if [ -e ./$filename ] ; then

     while read line ; do
          echo 'dev_pre= ' $line
          dev_pre=$line
     done < $filename

     # $dev_pre ( /= .or. == ) $dev     
     if [ $dev_pre -ne $dev ] ; then
          echo '$dev_pre -ne $dev'
          dev=$dev parameter_file=$parameter_file sh update_Parameters.sh
     else
          echo '$dev_pre -eq $dev'
     fi
else
     dev=$dev parameter_file=$parameter_file sh update_Parameters.sh
fi
#===========================================================================

#make clean
rm a.out

make

