rm *.gif
echo "======================================="
echo "            Start compile              "
echo "======================================="
mpiifort -mcmodel=large -g -o -traceback Objective_function.f90 EA22.f90 \
        EA22_component.f random.f90 \
        cmaes_subprogram.f90 cmaes.f90 mainprogram10_mpi.f90 \
        -mkl -o a.exe
#ifort Den_0927.f90 -mkl -mcmodel=largie
echo "======================================="
echo "        Start running program          "
echo "======================================="
mpirun -np 10 ./a.exe
./Create_Gif.sh	
