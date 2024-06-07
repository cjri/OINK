module load libs/gcc/11.4.0
module load libs/gsl/2.7/gcc-4.8.5
export PATH=/opt/gridware/depots/996bcebb/el7/pkg/compilers/gcc/11.4.0/bin/:$PATH
make -f Makefile PARALLEL=1
