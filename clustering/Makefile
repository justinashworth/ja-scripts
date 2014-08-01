CXX = g++

# compiler warning flags
WFLAGS = -Wall -W -Wextra -pedantic

# debug build flags
#DBFLAGS = -ggdb -pg
DBFLAGS = -ggdb
CXXFLAGS = $(WFLAGS) $(DBFLAGS)

# performance build flags
OFLAGS = -O3
CXXFLAGS = $(WFLAGS) $(OFLAGS)

EXE = fastcluster.bootstrap
OBJECTFILES = fastcluster.bootstrap.o

# external libraries
LDLIBS = -lstdc++
#LDFLAGS = -pg

# build targets
all: $(EXE)

$(EXE): $(OBJECTFILES)
#	$(CXX) $(OBJECTFILES) $(LDLIBS) -o $(EXE) $(LDFLAGS)
	$(CXX) $(OBJECTFILES) $(LDLIBS) -o $(EXE)

fastcluster.bootstrap.o: fastcluster.bootstrap.cpp
	$(CXX) fastcluster.bootstrap.cpp -c $(CXXFLAGS)

clean:
	-rm *.o $(EXE)

.PHONY: tags
tags:
	ctags *.cpp *.h
