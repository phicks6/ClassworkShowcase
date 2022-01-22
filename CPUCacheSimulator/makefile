all: memhier

clean:
	rm -f memhier
	rm -f trace.log
	rm -f *.txt
	

memhier: memhier.cpp
	g++ -std=c++11 -g -O0 -o memhier memhier.cpp

run:
	./memhier < trace.dat


