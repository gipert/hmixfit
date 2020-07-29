# Makefile
#
# Author: Luigi Pertoldi - pertoldi@pd.infn.it
# Created: Mon 6 May 2019

FLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++11 -g -O3 -I obj \
        $$(bat-config --cflags)
LIBS  = $$(bat-config --libs)
OBJ = obj/json.gch obj/utils.gch obj/GerdaFitter.o obj/gerda-fitter.o

all : dirs | gerda-fitter $(OBJ)

gerda-fitter : obj/GerdaFitter.o obj/gerda-fitter.o
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)

obj/gerda-fitter.o : gerda-fitter.cc
	$(CXX) $(FLAGS) -c -o $@ $< $(LIBS)

obj/GerdaFitter.o : GerdaFitter.cc GerdaFitter.hh
	$(CXX) $(FLAGS) -c -o $@ $< $(LIBS)

obj/%.gch : %.hpp
	$(CXX) $(FLAGS) -c -o $@ $^ $(LIBS)

dirs :
	mkdir -p obj

.PHONY : install clean dirs

clean :
	-rm -rf obj gerda-fitter

install : gerda-fitter
	install -d $(PREFIX)/bin
	install $< $(PREFIX)/bin
