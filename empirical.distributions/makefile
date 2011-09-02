CXX = g++

# compiler warning flags
WFLAGS = -Wall -W -Wextra -pedantic

# debug build flags
DBFLAGS = -ggdb
CXXFLAGS = $(WFLAGS) $(DBFLAGS)

# performance build flags
OFLAGS = -O3
#CXXFLAGS = $(WFLAGS) $(OFLAGS)

EXE = emp.distrib
OBJECTFILES = emp.distrib.o

# external libraries
LDLIBS = -lstdc++

# build targets
all: $(EXE)
$(EXE): $(OBJECTFILES)
	$(CXX) $(OBJECTFILES) $(LDLIBS) -o $(EXE)

clean:
	-rm *.o $(EXE)

.PHONY: tags
tags:
	ctags *.cpp *.h
