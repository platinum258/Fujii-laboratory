
# Parameters
Number_Loop_Frequency=0
Number_Loop_IncidentAngle=0
Number_Loop_DielectricConstant=0
Number_Loop_Position_Source_X=0
Number_Loop_Position_Source_Y=0

Flag_Plot_Configuration=1
Flag_Plot_Mesh_Configuration=0
Flag_Output_ElectricField=0
Flag_Plot_EFD_Configuration=1
Flag_Plot_EFD_OC_PoyntingVector=1

Flag_Polarization=0

if [ $dev -eq 0 ] ; then

     File_Name_Optimal=FEM_Data_Interval_001000
     File_Name_Normalization=FEM_Data_Interval_001000
     File_Name_Reference=FEM_Data_Interval_001000
     Pass_InputData_Optimal='../h_2xt'
     Pass_InputData_Normalization='../h_2xt'
     Pass_InputData_Reference='../h_2xt'

     Flag_Theoretical_Solution=0
     Flag_Reference_Device=1


elif [ $dev -eq 3 ] ; then

     File_Name_Optimal=Normalization
     File_Name_Normalization=Normalization
     File_Name_Reference=Reference
     Pass_InputData_Optimal='../FEM_Data/'
     Pass_InputData_Normalization='../FEM_Data/'
     Pass_InputData_Reference='../FEM_Data/'

     Flag_Theoretical_Solution=0
     Flag_Reference_Device=1

elif [ $dev -eq 10 ] ; then

     File_Name_Optimal=FEM_Data_Interval_00002000250
     File_Name_Normalization=FEM_Data_Interval_00002000250
     File_Name_Reference=FEM_Data_Interval_00002000250
     Pass_InputData_Optimal='../../../Optimization/ver1_max/result/'
     Pass_InputData_Normalization='../../../Optimization/ver1_max/result/'
     Pass_InputData_Reference='../../../Optimization/ver1_max/result/'

     Flag_Theoretical_Solution=0
     Flag_Reference_Device=1

elif [ $dev -eq 20 ] ; then
     File_Name_Optimal=FEM_Data_xmean_Convergence_004958
     File_Name_Normalization=FEM_Data_xmean_Convergence_004958
     File_Name_Reference=FEM_Data_xmean_Convergence_004958
     Pass_InputData_Optimal='../../../Optimization/tau_1e-2/result/'
     Pass_InputData_Normalization='../../../Optimization/tau_1e-2/result/'
     Pass_InputData_Reference='../../../Optimization/tau_1e-2/result/'

     Flag_Theoretical_Solution=0
     Flag_Reference_Device=1

elif [ $dev -eq 21 ] ; then
     File_Name_Optimal=FEM_Data_xmean_Convergence_004958
     File_Name_Normalization=FEM_Data_xmean_Convergence_004958
     File_Name_Reference=FEM_Data_xmean_Convergence_004958
     Pass_InputData_Optimal='../../../Optimization/tau_1e-2/result/'
     Pass_InputData_Normalization='../../../Optimization/tau_1e-2/result/'
     Pass_InputData_Reference='../../../Optimization/tau_1e-2/result/'

     Flag_Theoretical_Solution=0
     Flag_Reference_Device=1

elif [ $dev -eq 30 ] ; then
     File_Name_Optimal=FEM_Data_xmean_Convergence_003958
     File_Name_Normalization=${File_Name_Optimal}
     File_Name_Reference=${File_Name_Optimal}
     Pass_InputData_Optimal='../../../Optimization/tau_1e-2_restart_0.3_0.1/result/'
     Pass_InputData_Normalization=${Pass_InputData_Optimal}
     Pass_InputData_Reference=${Pass_InputData_Optimal}

     Flag_Theoretical_Solution=0
     Flag_Reference_Device=1

else
     echo 'dev=??? sh exe.sh'
     echo 'dev=0 : Thermal Cloak'
     echo 'dev=10 : Thermal Ground Cloak'

     echo 'setting.sh 50'
     exit 1
fi

File_Name_Result_Normalization=Result_Normalization.dat
File_Pass_Result_Normalization=../$Name_Directry_Normalization/
File_Pass_Reference_Field=../$Name_Directry_Reference/

sh clean_all.sh

cp -r ./Program ./$Name_Directry_Normalization
cp -r ./Program ./$Name_Directry_Optimal

rm ./Program/Parameters.f90
rm ./$Name_Directry_Normalization/Input_Data
rm ./$Name_Directry_Optimal/Input_Data


echo $Flag_Reference_Device

if [ $Flag_Reference_Device -eq 1 ] ; then
     cp -r ./Program ./$Name_Directry_Reference
     rm ./$Name_Directry_Reference/Input_Data

     # Reference
     echo $File_Name_Optimal >> ./$Name_Directry_Reference/Input_Data
     echo $File_Name_Normalization >> ./$Name_Directry_Reference/Input_Data
     echo $File_Name_Reference >> ./$Name_Directry_Reference/Input_Data
     echo "'"$Pass_InputData_Optimal"'" >> ./$Name_Directry_Reference/Input_Data
     echo "'"$Pass_InputData_Normalization"'" >> ./$Name_Directry_Reference/Input_Data
     echo "'"$Pass_InputData_Reference"'" >> ./$Name_Directry_Reference/Input_Data
     echo $Number_Loop_Frequency '# Number_Loop_Frequency' >> ./$Name_Directry_Reference/Input_Data
     echo $Number_Loop_IncidentAngle '# Number_Loop_IncidentAngle' >> ./$Name_Directry_Reference/Input_Data
     echo $Number_Loop_DielectricConstant '# Number_Loop_DielectricConstant' >> ./$Name_Directry_Reference/Input_Data
     echo $Number_Loop_Position_Source_X '# Number_Loop_Position_Source_X' >> ./$Name_Directry_Reference/Input_Data
     echo $Number_Loop_Position_Source_Y '# Number_Loop_Position_Source_Y' >> ./$Name_Directry_Reference/Input_Data
     echo '2 # Flag_Structure' >> ./$Name_Directry_Reference/Input_Data
     echo $Flag_Theoretical_Solution >> ./$Name_Directry_Reference/Input_Data
     echo 0 ' # Flag_Plot_Configuration' >> ./$Name_Directry_Reference/Input_Data
     #echo 1 ' # Flag_Plot_Configuration' >> ./$Name_Directry_Reference/Input_Data
     echo 0 ' # Flag_Plot_Mesh_Configuration' >> ./$Name_Directry_Reference/Input_Data
     #echo 1 '# Flag_Output_ElectricField' >> ./$Name_Directry_Reference/Input_Data
     echo 0 '# Flag_Output_ElectricField' >> ./$Name_Directry_Reference/Input_Data
     #echo 1 '# Flag_Plot_EFD_Configuration' >> ./$Name_Directry_Reference/Input_Data
     echo 0 '# Flag_Plot_EFD_Configuration' >> ./$Name_Directry_Reference/Input_Data
     #echo 1 '# Flag_Plot_EFD_OC_PoyntingVector' >> ./$Name_Directry_Reference/Input_Data
     echo 0 '# Flag_Plot_EFD_OC_PoyntingVector' >> ./$Name_Directry_Reference/Input_Data
     echo $File_Name_Result_Normalization '# File_Name_Result_Normalization' >> ./$Name_Directry_Reference/Input_Data
     echo "'"$File_Pass_Result_Normalization"'" '# File_Pass_Result_Normalization' >> ./$Name_Directry_Reference/Input_Data
     echo "'"$File_Pass_Reference_Field"'" '# File_Pass_Reference_Field' >> ./$Name_Directry_Reference/Input_Data
     echo $Flag_Reference_Device '# Flag_Reference_Device' >> ./$Name_Directry_Reference/Input_Data
fi

# Normalization
echo $File_Name_Optimal >> ./$Name_Directry_Normalization/Input_Data
echo $File_Name_Normalization >> ./$Name_Directry_Normalization/Input_Data
echo $File_Name_Reference >> ./$Name_Directry_Normalization/Input_Data
echo "'"$Pass_InputData_Optimal"'" >> ./$Name_Directry_Normalization/Input_Data
echo "'"$Pass_InputData_Normalization"'" >> ./$Name_Directry_Normalization/Input_Data
echo "'"$Pass_InputData_Reference"'" >> ./$Name_Directry_Normalization/Input_Data
echo $Number_Loop_Frequency '# Number_Loop_Frequency' >> ./$Name_Directry_Normalization/Input_Data
echo $Number_Loop_IncidentAngle '# Number_Loop_IncidentAngle' >> ./$Name_Directry_Normalization/Input_Data
echo $Number_Loop_DielectricConstant '# Number_Loop_DielectricConstant' >> ./$Name_Directry_Normalization/Input_Data
echo $Number_Loop_Position_Source_X '# Number_Loop_Position_Source_X' >> ./$Name_Directry_Normalization/Input_Data
echo $Number_Loop_Position_Source_Y '# Number_Loop_Position_Source_Y' >> ./$Name_Directry_Normalization/Input_Data
echo '1 # Flag_Structure' >> ./$Name_Directry_Normalization/Input_Data
echo $Flag_Theoretical_Solution >> ./$Name_Directry_Normalization/Input_Data
echo 0 ' # Flag_Plot_Configuration' >> ./$Name_Directry_Normalization/Input_Data
#echo 1 ' # Flag_Plot_Configuration' >> ./$Name_Directry_Normalization/Input_Data
echo 0 ' # Flag_Plot_Mesh_Configuration' >> ./$Name_Directry_Normalization/Input_Data
#echo 1 '# Flag_Output_ElectricField' >> ./$Name_Directry_Normalization/Input_Data
echo 0 '# Flag_Output_ElectricField' >> ./$Name_Directry_Normalization/Input_Data
#echo 1 '# Flag_Plot_EFD_Configuration' >> ./$Name_Directry_Normalization/Input_Data
echo 0 '# Flag_Plot_EFD_Configuration' >> ./$Name_Directry_Normalization/Input_Data
#echo 1 '# Flag_Plot_EFD_OC_PoyntingVector' >> ./$Name_Directry_Normalization/Input_Data
echo 0 '# Flag_Plot_EFD_OC_PoyntingVector' >> ./$Name_Directry_Normalization/Input_Data
echo $File_Name_Result_Normalization '# File_Name_Result_Normalization' >> ./$Name_Directry_Normalization/Input_Data
echo "'"$File_Pass_Result_Normalization"'" '# File_Pass_Result_Normalization' >> ./$Name_Directry_Normalization/Input_Data
echo "'"$File_Pass_Reference_Field"'" '# File_Pass_Reference_Field' >> ./$Name_Directry_Normalization/Input_Data
echo $Flag_Reference_Device '# Flag_Reference_Device' >> ./$Name_Directry_Normalization/Input_Data

# Optimal
echo $File_Name_Optimal >> ./$Name_Directry_Optimal/Input_Data
echo $File_Name_Normalization >> ./$Name_Directry_Optimal/Input_Data
echo $File_Name_Reference >> ./$Name_Directry_Optimal/Input_Data
echo "'"$Pass_InputData_Optimal"'" >> ./$Name_Directry_Optimal/Input_Data
echo "'"$Pass_InputData_Normalization"'" >> ./$Name_Directry_Optimal/Input_Data
echo "'"$Pass_InputData_Reference"'" >> ./$Name_Directry_Optimal/Input_Data
echo $Number_Loop_Frequency '# Number_Loop_Frequency' >> ./$Name_Directry_Optimal/Input_Data
echo $Number_Loop_IncidentAngle '# Number_Loop_IncidentAngle' >> ./$Name_Directry_Optimal/Input_Data
echo $Number_Loop_DielectricConstant '# Number_Loop_DielectricConstant' >> ./$Name_Directry_Optimal/Input_Data
echo $Number_Loop_Position_Source_X '# Number_Loop_Position_Source_X' >> ./$Name_Directry_Optimal/Input_Data
echo $Number_Loop_Position_Source_Y '# Number_Loop_Position_Source_Y' >> ./$Name_Directry_Optimal/Input_Data
echo '0 # Flag_Structure' >> ./$Name_Directry_Optimal/Input_Data
echo $Flag_Theoretical_Solution >> ./$Name_Directry_Optimal/Input_Data
echo $Flag_Plot_Configuration ' # Flag_Plot_Configuration' >> ./$Name_Directry_Optimal/Input_Data
echo $Flag_Plot_Mesh_Configuration ' # Flag_Plot_Mesh_Configuration' >> ./$Name_Directry_Optimal/Input_Data
echo $Flag_Output_ElectricField '# Flag_Output_ElectricField' >> ./$Name_Directry_Optimal/Input_Data
echo $Flag_Plot_EFD_Configuration '# Flag_Plot_EFD_Configuration' >> ./$Name_Directry_Optimal/Input_Data
echo $Flag_Plot_EFD_OC_PoyntingVector '# Flag_Plot_EFD_OC_PoyntingVector' >> ./$Name_Directry_Optimal/Input_Data
echo $File_Name_Result_Normalization '# File_Name_Result_Normalization' >> ./$Name_Directry_Optimal/Input_Data
echo "'"$File_Pass_Result_Normalization"'" '# File_Pass_Result_Normalization' >> ./$Name_Directry_Optimal/Input_Data
echo "'"$File_Pass_Reference_Field"'" '# File_Pass_Reference_Field' >> ./$Name_Directry_Optimal/Input_Data
echo $Flag_Reference_Device '# Flag_Reference_Device' >> ./$Name_Directry_Optimal/Input_Data 


