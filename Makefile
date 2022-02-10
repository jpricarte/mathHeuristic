GUROBI_DIR=opt/gurobi950/linux64
FLAGS=-std=c++11 -Wall -Wextra -lemon -fsplit-stack
INCLUDE=-I/$(GUROBI_DIR)/include/
LIBRARY=-L/$(GUROBI_DIR)/lib/ -lgurobi_c++ -lgurobi95

heuristic: heuristic.cpp
	g++ $(FLAGS) heuristic.cpp -o heuristic $(INCLUDE) $(LIBRARY) 

rpheuro: rpheuro.cpp
	g++ $(FLAGS) rpheuro.cpp -o rpheuro $(INCLUDE) $(LIBRARY)

heuristic2: heuristicBAK.cpp
	g++ $(FLAGS) heuristicBAK.cpp -o heuristic $(INCLUDE) $(LIBRARY) 