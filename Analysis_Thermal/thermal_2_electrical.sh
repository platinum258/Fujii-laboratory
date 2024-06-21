grep -l 'Thermal' ./Program/*.f90 | xargs sed -i.bak -e 's/Thermal/Electrical/g'
grep -l 'Thermal' ./Program/Parameters/*.f90 | xargs sed -i.bak -e 's/Thermal/Electrical/g'
grep -l 'Thermal' ./Program/makefile | xargs sed -i.bak -e 's/Thermal/Electrical/g'
grep -l 'Thermal' ./exe_all.sh | xargs sed -i.bak -e 's/Thermal/Electrical/g'
grep -l 'Thermal' ./all_compile.sh | xargs sed -i.bak -e 's/Thermal/Electrical/g'
rename Thermal Electrical ./Program/*Thermal*.f90
rename Thermal Electrical ./Program/Parameters/*Thermal*.f90
#ls ./Program/Parameters/*Electrical*.f90
grep -l 'Heat' ./Program/*.f90 | xargs sed -i.bak -e 's/Heat/Electric/g'
grep -l 'Heat' ./Program/Parameters/*.f90 | xargs sed -i.bak -e 's/Heat/Electric/g'
grep -l 'Heat' ./Program/makefile | xargs sed -i.bak -e 's/Heat/Electric/g'
grep -l 'Heat' ./exe_all.sh | xargs sed -i.bak -e 's/Heat/Electric/g'
grep -l 'Heat' ./all_compile.sh | xargs sed -i.bak -e 's/Heat/Electric/g'
grep -l 'Electric_Flux' ./Program/*.f90 | xargs sed -i.bak -e 's/Electric_Flux/Direct_Current/g'
rename Heat Electric ./Program/*Heat*.f90
rename Heat Electric ./Program/Parameters/*Heat*.f90
#ls ./Program/Parameters/*Electric*.f90
