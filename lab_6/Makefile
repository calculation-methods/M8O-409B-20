matrix.o: math/matrix.cpp math/matrix.hpp 
	g++ -c math/matrix.cpp -o tmp/matrix.o

lab.o: lab6.cpp lab6.hpp 
	g++ -c lab6.cpp -o tmp/lab.o

solver.o: math/solver.cpp math/solver.hpp
	g++ -c math/solver.cpp -o tmp/solver.o

main.o: main.cpp
	g++ -c main.cpp -o tmp/main.o

prog: main.o matrix.o solver.o lab.o
	g++ tmp/matrix.o tmp/main.o tmp/solver.o tmp/lab.o -o prog
