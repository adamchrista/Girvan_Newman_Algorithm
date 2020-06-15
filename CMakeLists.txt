cmake_minimum_required(VERSION 3.12)
project(20s-pa02-doover-adamchrista)
set (EXE_NAME 20s_pa02_doover_adamchrista)

##########################################
# Students, edit the items below.  In the quotes,
# put the names of the input files as they should appear
# in the argv param to main.
# Example:
#   set (input01 "train_data.csv")
#   set (input02 "train_target.csv")
# etc....

set(input01 "cmd.txt")
set(input02 "graph.in")
set(input03 "graph.txt")


#
##########################################

# Copy data files to build dir
foreach(file IN LISTS input01 input02 input03)
    configure_file(${file} ${file} COPYONLY)
endforeach()

set(CMAKE_CXX_STANDARD 17)
add_executable(${EXE_NAME} main.cpp Edge.h Graph.h Node.h graph.in cmd.txt graph.txt)
