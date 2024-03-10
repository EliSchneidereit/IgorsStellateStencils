igor.exe : igor.cpp
	g++ -O3 $< -o $@ -lboost_program_options

all : igor.exe
