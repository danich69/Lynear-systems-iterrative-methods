Jacobi: prog
	./prog Jacobi
Zeidel: prog
	./prog Zeidel
Relax: prog
	./prog Relax
prog: main.o Linear_system_iter.o
	ifort $^ -o $@ -debug
main.o: main.f90 linear_system_iter.mod
	ifort $^ -c
linear_system_iter.mod Linear_system_iter.o: Linear_system_iter.f90
	ifort Linear_system_iter.f90 -c
clean: 
	rm -f *.o *mod cr prog
