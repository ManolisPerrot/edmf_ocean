FC = gfortran
CC = gcc
FFLAGS = -O4 -fdefault-real-8 -fdefault-double-8
f2pyFFLAGS = "-O4 -fdefault-real-8 -fdefault-double-8"
#FFLAGS = -g -O0 -fcheck=bounds -Wall -Wextra -fbacktrace -finit-real=snan -finit-integer=8888 -finit-character=90
#f2pyFFLAGS = "-g -O0 -fcheck=bounds -Wall -Wextra -fbacktrace -finit-real=snan -finit-integer=8888 -finit-character=90"

.PHONY: doc clean clean_all

GENERATED_FILES = scm_oce_new.f90, scm_mfc_new.f90, scm_keps_new.f90, scm_tke_new.f90

all : scm_par.o scm_oce_new.so scm_tke_new.so scm_mfc_new.so scm_keps_new.so 
	@cp -f scmoce_new.cpython-311-x86_64-linux-gnu.so ../F2PY/scmoce_new.so 
	@cp -f scmtke_new.cpython-311-x86_64-linux-gnu.so ../F2PY/scmtke_new.so
	@cp -f scmmfc_new.cpython-311-x86_64-linux-gnu.so ../F2PY/scmmfc_new.so
	@cp -f scmkeps_new.cpython-311-x86_64-linux-gnu.so ../F2PY/scmkeps_new.so
	@/bin/rm -f *.o *.mod
	@echo
	@echo ===================================================
	@echo       Fortran SCM tools OK
	@echo ===================================================
	@echo

scm_oce_new.so : scm_par.o scm_oce_new.f90
	f2py -c scm_oce_new.f90 scm_par.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scmoce_new

scm_tke_new.so : scm_par.o scm_tke_new.f90
	f2py -c scm_tke_new.f90 scm_par.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scmtke_new

scm_mfc_new.so : scm_par.o scm_mfc_new.f90
	f2py -c scm_mfc_new.f90 scm_par.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scmmfc_new

scm_keps_new.so : scm_par.o scm_keps_new.f90
	f2py -c scm_keps_new.f90 scm_par.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scmkeps_new

scm_par.o : scm_par.f90
	$(FC) $(FFLAGS) -o scm_par.o -c scm_par.f90

clean:
	/bin/rm -f *.o *.mod ../F2PY/scmkeps_new.so  ../F2PY/scmmfc_new.so  ../F2PY/scmoce_new.so  ../F2PY/scmtke_new.so *.pyc *.msg* *~
clean_all:
	/bin/rm -f *.o *.mod ../F2PY/scmkeps_new.so  ../F2PY/scmmfc_new.so  ../F2PY/scmoce_new.so  ../F2PY/scmtke_new.so *.pyc *.msg* *~
