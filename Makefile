CC	      = gcc
#CC_FLAGS        = -g -O2 -Wall -D_GLIBCXX_DEBUG -pg
#LD_FLAGS	= -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl -lgslcblas -lm -lstdc++ -pg -g
CC_FLAGS        = -g -O3 -Wall 
LD_FLAGS	= -lgsl -lgslcblas -lm -lstdc++  -g

BAS		= basicmodel.o io.o utilities.o 

basic: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o oink  $(LD_FLAGS)
basicmodel.o: basicmodel.cpp
	$(CC) $(CC_FLAGS) -c basicmodel.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp
pseudorandom.o: pseudorandom.cpp
	 $(CC) $(CC_FLAGS) -c pseudorandom.cpp

