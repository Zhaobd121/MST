CXX = g++
CXXFLAGS = -std=c++17 -Wall -O3

# Target executable
TARGET = program

# Source files
SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)

# Header files
HEADERS = graph.h hash.hpp disjointSet.h heap.hpp timer.h utility.hpp FCI.h tree.h tree1.h

# Build target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)
