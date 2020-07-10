default: crossgrad_solver crossgrad.so

crossgrad.so : crossgrad_mod.f90
	f2py3 --overwrite-signature $< -m crossgrad -h crossgrad.pyf
	f2py3 -c $< crossgrad.pyf

crossgrad_solver : crossgrad_mod.o crossgrad_solver.o
	gfortran -o $@ $^

%.o : %.f90
	gfortran -c -Wall -O2 $<

clean:
	rm -f *~
	rm -f *.pyf *.mod *.o *.so
