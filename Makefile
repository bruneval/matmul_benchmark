FCFLAGS= -fpp
SCALAPACK_LDFLAGS=${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl -liomp5

#SCALAPACK=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm

SCALAPACK= -DSCALAPACK ${SCALAPACK_LDFLAGS} -mkl=parallel

all: matmul_pdgemm matmul_dgemm matmul_sgemm

matmul_pdgemm: matmul_pdgemm.f90
	mpifort ${FCFLAGS} -o matmul_pdgemm matmul_pdgemm.f90  -O2 -xhost ${SCALAPACK}

matmul_dgemm: matmul_pdgemm.f90
	ifort ${FCFLAGS} -o matmul_dgemm matmul_pdgemm.f90 -O2 -xhost -mkl=parallel 

matmul_sgemm: matmul_pdgemm.f90
	ifort ${FCFLAGS} -o matmul_sgemm matmul_pdgemm.f90 -O2 -xhost -mkl=parallel -DSINGLE_PREC

clean: 
	rm -f matmul_pdgemm matmul_dgemm matmul_sgemm
