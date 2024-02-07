CC	      = gcc
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS	= -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl -lgslcblas -lm -lstdc++ 
BAS		= basicmodel.o io.o utilities.o pseudorandom.o

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

