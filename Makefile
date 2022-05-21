GUROBI_DIR=opt/gurobi951/linux64
FLAGS=-std=c++11 -Wall -Wextra -lemon -fsplit-stack
INCLUDE=-I/$(GUROBI_DIR)/include/
LIBRARY=-L/$(GUROBI_DIR)/lib/ -lgurobi_c++ -lgurobi95

heuristic: heuristic.cpp
	g++ $(FLAGS) heuristic.cpp -o heuristic $(INCLUDE) $(LIBRARY) 

lp: lp.cpp
	g++ $(FLAGS) lp.cpp -o lp $(INCLUDE) $(LIBRARY)