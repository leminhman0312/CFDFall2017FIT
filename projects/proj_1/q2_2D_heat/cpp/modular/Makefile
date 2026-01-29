CXX = g++
CXXFLAGS = -O2 -std=c++17 -Iinclude

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)

a.out: $(OBJ)
	$(CXX) $(OBJ) -o a.out

clean:
	rm -f src/*.o a.out
