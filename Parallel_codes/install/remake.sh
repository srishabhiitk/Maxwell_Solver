rm -rf build
rm -f ../aether
mkdir build
cd build
cmake ../
make
cd ../../
#./aether
mpirun -np 8 ./maxwell
