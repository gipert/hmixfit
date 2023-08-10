# Makefile
#
# Author: Luigi Pertoldi - pertoldi@pd.infn.it
# Created: Mon 6 May 2019

FLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++11 -g -O3 -I build \
        $$(bat-config --cflags)
LIBS  = $$(bat-config --libs)
OBJ = build/HMixFit.o build/hmixfit.o build/json.gch build/utils.gch

all : dirs | build/hmixfit $(OBJ)

build/hmixfit : $(OBJ)
	$(CXX) $(FLAGS) -o $@ build/HMixFit.o build/hmixfit.o $(LIBS)

build/hmixfit.o : src/hmixfit.cc src/json.hpp
	$(CXX) $(FLAGS) -c -o $@ $< $(LIBS)

build/HMixFit.o : src/HMixFit.cc src/HMixFit.hh src/utils.hpp src/json.hpp
	$(CXX) $(FLAGS) -c -o $@ $< $(LIBS)

build/%.gch : src/%.hpp
	$(CXX) $(FLAGS) -c -o $@ $^ $(LIBS)

dirs :
	mkdir -p build

.PHONY : install clean dirs

clean :
	-rm -rf build

install : dirs | build/hmixfit
	install -d $(PREFIX)/bin
	install build/hmixfit $(PREFIX)/bin
