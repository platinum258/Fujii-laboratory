
rm Device_Number
echo $dev >> Device_Number
rm Parameters.o
rm Parameters.f90
cp ./Parameters/$parameter_file ./Parameters.f90
echo '    '
echo '    ================================================= '
echo "     $parameter_file  --> Parameters.f90"
echo '    ================================================= '
echo '    '

