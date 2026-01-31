CXX := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -Iinclude

TARGET := 2D_heat

SRC := $(wildcard src/*.cpp)
OBJ := $(patsubst src/%.cpp, build/%.o, $(SRC))

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p build

clean:
	rm -f $(TARGET) build/*.o data/*.txt data/*.dat plot/*.png

run: $(TARGET)
	./$(TARGET)
