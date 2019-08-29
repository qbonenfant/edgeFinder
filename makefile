# Example commands:
#   make (build in release mode)
#   make debug (build in debug mode)
#   make clean (deletes *.o files)
#   make distclean (deletes *.o files and the *.so files)
#   make CXX=g++-5 (build with a particular compiler)
#   make CXXFLAGS="-Werror -g3" (build with particular compiler flags)


# CXX and CXXFLAGS can be overridden by the user.
CXX         ?= g++-5
CXXFLAGS    ?= -Wall -Wextra -pedantic -mtune=native 

# These flags are required for the build to work.
FLAGS        = -std=c++14

# Different debug/optimisation levels for debug/release builds.
DEBUGFLAGS   = -DSEQAN_ENABLE_DEBUG=1 -g
RELEASEFLAGS = -O3 -DNDEBUG -w

TARGET       = ./edgeFinder
SHELL        = /bin/sh
SOURCES      = $(shell find ./ -name "*.cpp")
HEADERS      = $(shell find ./ -name "*.h")
OBJECTS      = $(SOURCES:.cpp=.o)

# Linux needs '-soname' while Mac needs '-install_name'
PLATFORM     = $(shell uname)
ifeq ($(PLATFORM), Darwin)
SONAME       = -install_name
OMP 		 = 
LRT 		 = 

else
SONAME       = -soname
OMP 		 = -fopenmp
LRT 		 = -lrt
endif


all: $(TARGET)

.PHONY: release
release: FLAGS+=$(RELEASEFLAGS)
release: $(TARGET)

.PHONY: debug
debug: FLAGS+=$(DEBUGFLAGS)
debug: $(TARGET)


$(TARGET): $(OBJECTS)
	$(CXX) $(FLAGS) $(CXXFLAGS)  -Wl,$(SONAME),$(TARGET) -o $(TARGET) $(OBJECTS) $(LRT) -l pthread

clean:
	$(RM) $(OBJECTS)

distclean: clean
	$(RM) $(TARGET) 

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) $(CXXFLAGS) -c -o $@ $<
