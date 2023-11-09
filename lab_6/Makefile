matrix.o: math/matrix.cpp math/matrix.hpp
	g++ -c math/matrix.cpp -o tmp/matrix.o

solver.o: math/solver.cpp math/solver.hpp
	g++ -c math/solver.cpp -o tmp/solver.o

lab_6.o: lab_6.cpp lab_6.hpp
	g++ -c lab_6.cpp -o tmp/lab_6.o

main.o: main.cpp
	g++ -c main.cpp -o tmp/main.o

prog: main.o matrix.o solver.o
	g++ tmp/matrix.o tmp/main.o tmp/solver.o tmp/lab_6.o -o prog
