BUILD_DIR	?= build/
SRC_DIR		?= src/


CC		= g++
CC_FLAGS        = -g3 -O3 -Wall -std=c++11 -pthread  -g $(if $(PARALLEL), -fopenmp) -I /opt/homebrew/Cellar/gsl/2.7.1/include/ $(if $(GSL_INCLUDE), -I $(GSL_INCLUDE))
LD_FLAGS	= -L/opt/homebrew/Cellar/gsl/2.7.1/lib $(if $(GSL_LIB), -L $(GSL_LIB)) -lgsl -lgslcblas -lm -g 
LD_FLAGS_TEST   = $(if $(GSL_LIB), -L $(GSL_LIB)) -lgtest -lgtest_main -pthread -lgsl -lgslcblas -g
BAS		= $(BUILD_DIR)/basicmodel.o $(BUILD_DIR)/io.o $(BUILD_DIR)/utilities.o 
TEST            = $(BUILD_DIR)/Test.o $(BUILD_DIR)/utilities.o $(BUILD_DIR)/io.o

basic: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o oink  $(LD_FLAGS)
test: $(TEST) 
	$(CC) $(CC_FLAGS) $(TEST) -o $(SRC_DIR)/test $(LD_FLAGS_TEST)
$(BUILD_DIR)/Test.o: ${SRC_DIR}/Test.cpp | $(BUILD_DIR)
	$(CC) $(CC_FLAGS) -c $(SRC_DIR)/Test.cpp -o $(BUILD_DIR)/Test.o
$(BUILD_DIR)/basicmodel.o: $(SRC_DIR)/basicmodel.cpp | $(BUILD_DIR)
	$(CC) $(CC_FLAGS) -c $(SRC_DIR)/basicmodel.cpp -o $(BUILD_DIR)/basicmodel.o  
$(BUILD_DIR)/io.o: $(SRC_DIR)/io.cpp | $(BUILD_DIR)
	$(CC) $(CC_FLAGS) -c $(SRC_DIR)/io.cpp -o $(BUILD_DIR)/io.o
$(BUILD_DIR)/utilities.o: $(SRC_DIR)/utilities.cpp | $(BUILD_DIR)
	$(CC) $(CC_FLAGS) -c $(SRC_DIR)/utilities.cpp -o $(BUILD_DIR)/utilities.o
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

