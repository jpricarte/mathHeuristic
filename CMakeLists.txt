cmake_minimum_required(VERSION 3.23)
project(math_heuristic)

set(CMAKE_CXX_STANDARD 17)

# declara header gurobi directory
set(GUROBI_INCLUDE_DIRECTORY "/opt/gurobi951/linux64/include")

# adicionar header gurobi
include_directories(${GUROBI_INCLUDE_DIRECTORY})

# declara library gurobi directory
set(GUROBI_LIB_DIRECTORY "/opt/gurobi951/linux64/lib")

# adicionar library gurobi
link_directories(${GUROBI_LIB_DIRECTORY})

add_executable(math_heuristic main.cpp ${GUROBI_INCLUDE_DIRECTORY}/gurobi_c++.h instance.h heuristic.h auxiliar.h)
target_link_libraries(math_heuristic emon)
target_link_libraries(math_heuristic gurobi95)
target_link_libraries(math_heuristic gurobi_c++)
