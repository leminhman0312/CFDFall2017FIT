CXX := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra

TARGET := 1D_heat

SRC := main.cpp helpers.cpp methods.cpp driver.cpp
OBJ := $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

cleanall: clean
	rm -rf plot compare 

run: $(TARGET)
	./$(TARGET)
