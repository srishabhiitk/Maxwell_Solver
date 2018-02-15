rm -rf build
rm -f ../maxwell
mkdir build
cd build
SCOREP_WRAPPER_OFF=true cmake ../
make  SCOREP_WRAPPER_ARGS="--mpp=mpi --verbose" \
      SCOREP_WRAPPER_FLAGS="-O3"
cd ../../
#./aether
mpirun -np 4 ./maxwell
