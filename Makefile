igor.exe : igor.cpp
	g++ -O3 $< -o $@

all : igor.exe
