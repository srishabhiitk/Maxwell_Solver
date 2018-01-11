rm -f ../aether
cd build
make
cd ../../
#./aether
#mpirun -np 16 ./aether
mpirun -np 8 ./maxwell
#scp output/Output_Data.dat sam:
