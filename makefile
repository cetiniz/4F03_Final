CC=mpiCC
CFLAGS=-c -g -Wall -std=c++11 -fopenmp
CXXFLAGS=-fopenmp
LDFLAGS=-lm
OBJECTS=$(SOURCES:.cpp=.o)

SOURCES=main.cpp savebmp.c

EXECUTABLE=x.project

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $(CXXFLAGS) $< -o $@

clean:
	rm *.o x.*
