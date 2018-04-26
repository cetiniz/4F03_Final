CC=mpiCC
CFLAGS=-c -g -Wall -std=c++11 -fopenmp
CXXFLAGS=-fopenmp
LDFLAGS=-lm
OBJECTS=$(SOURCES:.cpp=.o)

SOURCES=main.cpp savebmp.c

EXECUTABLE=project.x

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o *.x
