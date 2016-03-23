# makefile ffor GST project

# external tool locations
BOOST_HOME  = /usr/local/boost_1_60_0
GUROBI_HOME = /opt/gurobi650/linux64

SOURCE   = src
INCL     = include
INCLUDE  = -I./include -I$(BOOST_HOME) -I$(GUROBI_HOME)/include
BUILD    = build


OBJ_GCC = $(BUILD)/DataTypes.o $(BUILD)/GST.o $(BUILD)/main.o 


CC     :=  g++
CFLAGS  = -std=c++11 -O3 
OFLAGS := -g3
OBJ     = $(OBJ_GCC)
LIBS    = -L$(GUROBI_HOME)/lib/ -L$(BOOST_HOME)/stage/lib -lgurobi_c++ -lgurobi65 -lm -lboost_filesystem -lboost_system


all: GST_solver
	@echo  Project GST_solver has been compiled.

GST_solver: $(OBJ)
	@$(CC) $(CFLAGS) $(OBJ) -o GST_solver $(LIBS)

$(BUILD)/%.o : $(SOURCE)/%.cpp
	@$(CC) $(CFLAGS) $(OFLAGS) $(INCLUDE) $< -c -o $@ 



test:
	@echo $(CC)



clean:
	@rm -f $(BUILD)/*.o GST_solver


	
	
	
	
