rm *.gif
echo "======================================="
echo "            Start compile              "
echo "======================================="
mpif90 -mcmodel=large -g Objective_function.f90 EA22.f90 \
        EA22_component.f random.f90 \
        cmaes_subprogram.f90 cmaes.f90 \
        analyze_thermal_distribution.f90 Parameters.f90\
        Format_Global_aa_DMUMPS.f90 Format_Local_Global_DMUMPS.f90 \
        Implement_Dirichlet_BC_real_CSRF.f90 Implement_Neumann_BC_complex_CSRF.f90 \
        MA57.f90 MA57_Component.f \
        Output_Error.f90 mainprogram10_mpi.f90 \
        -mkl -o a.exe
#ifort Den_0927.f90 -mkl -mcmodel=largie
echo "======================================="
echo "        Start running program          "
echo "======================================="
mpirun -np 10 ./a.exe
./Create_Gif.sh	
