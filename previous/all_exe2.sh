
echo "======================================="
echo "            Start compile              "
echo "======================================="
ifort -mcmodel=medium Objective_function.f90 EA22.f90 \
        EA22_component.f random.f90 \
        femprogram4_1.f90 \
        -mkl -o a.exe
#ifort Den_0927.f90 -mkl -mcmodel=largie
echo "======================================="
echo "        Start running program          "
echo "======================================="
./a.exe

