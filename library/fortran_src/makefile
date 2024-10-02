FC = gfortran
CC = gcc
FFLAGS = -O4 -fdefault-real-8 -fdefault-double-8
f2pyFFLAGS = "-O4 -fdefault-real-8 -fdefault-double-8"
#FFLAGS = -g -O0 -fcheck=bounds -Wall -Wextra -fbacktrace -finit-real=snan -finit-integer=8888 -finit-character=90
#f2pyFFLAGS = "-g -O0 -fcheck=bounds -Wall -Wextra -fbacktrace -finit-real=snan -finit-integer=8888 -finit-character=90"

.PHONY: doc clean clean_all

GENERATED_FILES = scm_oce.f90, scm_tke.f90, scm_mfc.f90, scm_keps.f90

all : scm_par.o scm_oce.so scm_tke.so scm_mfc.so scm_keps.so
	@cp -f scmoce.cpython-311-x86_64-linux-gnu.so ../F2PY/scmoce.so 
	@cp -f scmtke.cpython-311-x86_64-linux-gnu.so ../F2PY/scmtke.so
	@cp -f scmmfc.cpython-311-x86_64-linux-gnu.so ../F2PY/scmmfc.so
	@cp -f scmkeps.cpython-311-x86_64-linux-gnu.so ../F2PY/scmkeps.so
	@/bin/rm -f *.o *.mod
	@echo
	@echo ===================================================
	@echo       Fortran SCM tools OK
	@echo ===================================================
	@echo

scm_oce.so : scm_par.o scm_oce.f90
	f2py -c scm_oce.f90 scm_par.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scmoce

scm_tke.so : scm_par.o scm_tke.f90
	f2py -c scm_tke.f90 scm_par.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scmtke

scm_mfc.so : scm_par.o scm_mfc.f90
	f2py -c scm_mfc.f90 scm_par.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scmmfc

scm_keps.so : scm_par.o scm_keps.f90
	f2py -c scm_keps.f90 scm_par.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scmkeps

scm_par.o : scm_par.f90
	$(FC) $(FFLAGS) -o scm_par.o -c scm_par.f90

clean:
	/bin/rm -f *.o *.mod ../F2PY/*.so *.pyc *.msg* *~
clean_all:
	/bin/rm -f *.o *.mod *.so *.pyc *.msg* *~
