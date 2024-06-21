
sh clean_all.sh

if [ -z $dev ] ; then
     shell_name=exe.sh line_number=3 sh dev_echo.sh 
     exit 1
fi

if [ $dev -eq 0 ] ; then
     parameter_file='Parameters_Thermal_Cloak.f90'

     Name_Directry_Optimal=Optimal_configuration
     Name_Directry_Normalization=Normalization
     Name_Directry_Reference=Reference

elif [ $dev -eq 3 ] ; then
     parameter_file='Parameters_Thermal_Cloak.f90'

     Name_Directry_Optimal=Optimal_configuration
     Name_Directry_Normalization=Normalization
     Name_Directry_Reference=Reference

elif [ $dev -eq 10 ] ; then
     parameter_file='Parameters_Thermal_Ground_Cloak.f90'

     Name_Directry_Optimal=Optimal_configuration
     Name_Directry_Normalization=Normalization
     Name_Directry_Reference=Reference

elif [ $dev -eq 20 ] ; then
     parameter_file='Parameters_Thermal_Cloak_Concentrator.f90'

     Name_Directry_Optimal=Optimal_configuration
     Name_Directry_Normalization=Normalization
     Name_Directry_Reference=Reference

elif [ $dev -eq 21 ] ; then
     parameter_file='Parameters_Thermal_Cloak_Concentrator_Source.f90'

     Name_Directry_Optimal=Optimal_configuration
     Name_Directry_Normalization=Normalization
     Name_Directry_Reference=Reference

elif [ $dev -eq 30 ] ; then
     parameter_file='Parameters_Thermal_Reversal.f90'

     Name_Directry_Optimal=Optimal_configuration
     Name_Directry_Normalization=Normalization
     Name_Directry_Reference=Reference

else

     shell_name=exe.sh line_number=67 sh dev_echo.sh 
     exit 1

fi

dev=$dev \
Name_Directry_Optimal=$Name_Directry_Optimal \
Name_Directry_Normalization=$Name_Directry_Normalization \
Name_Directry_Reference=$Name_Directry_Reference \
sh setting.sh

#===
if [ $dev -eq 0 -o $dev -eq 3 -o $dev -eq 10 -o $dev -eq 20 -o $dev -eq 21 -o $dev -eq 30 ] ; then
   cd $Name_Directry_Reference
   sh clean.sh
   cd ../
fi

cd ./$Name_Directry_Normalization
sh clean.sh

cd ../$Name_Directry_Optimal
sh clean.sh

cd ..
#===
if [ $dev -eq 0 -o $dev -eq 3 -o $dev -eq 10 -o $dev -eq 20 -o $dev -eq 21 -o $dev -eq 30 ] ; then
   cd ./$Name_Directry_Reference
   dev=$dev parameter_file=$parameter_file sh exe.sh
   cd ..
fi

cd ./$Name_Directry_Normalization
dev=$dev parameter_file=$parameter_file sh exe.sh
cd ..

cd ./$Name_Directry_Optimal
dev=$dev parameter_file=$parameter_file sh exe.sh
cd ..

