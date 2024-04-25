echo "======================================="
echo "            Start compile              "
echo "======================================="
ifort -mcmodel=medium 23W4038G_femprogram.f90 \
        -mkl -o a.exe
#ifort Den_0927.f90 -mkl -mcmodel=largie
echo "======================================="
echo "        Start running program          "
echo "======================================="
./a.exe