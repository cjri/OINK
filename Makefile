CC	      = g++
CC_FLAGS        = -g3 -O3 -Wall -pthread  -g $(if $(PARALLEL),-fopenmp)
LD_FLAGS	= -lgsl -lgslcblas -lm -lstdc++  -g 
LD_FLAGS_TEST   = -lgtest -lgtest_main -pthread -lgsl -lgslcblas -g
BAS		= basicmodel.o io.o utilities.o 
TEST            = Test.o utilities.o io.o
TEST_RNG        = test_rng.o utilities.o io.o

basic: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o oink  $(LD_FLAGS)
test: $(TEST)
	$(CC) $(CC_FLAGS) $(TEST) -o test $(LD_FLAGS_TEST)
test_rng: $(TEST_RNG)
	$(CC) $(CC_FLAGS) $(TEST_RNG) -o test $(LD_FLAGS)
test_rng.o: test_rng.cpp
	$(CC) $(CC_FLAGS) -c test_rng.cpp
Test.o: Test.cpp
	$(CC) $(CC_FLAGS) -c Test.cpp
basicmodel.o: basicmodel.cpp
	$(CC) $(CC_FLAGS) -c basicmodel.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp
pseudorandom.o: pseudorandom.cpp
	 $(CC) $(CC_FLAGS) -c pseudorandom.cpp

