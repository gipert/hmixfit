# Makefile
#
# Author: Luigi Pertoldi - pertoldi@pd.infn.it
# Created: Mon 6 May 2019

FLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++11 -g -O3 -I obj \
        $$(bat-config --cflags)
LIBS  = $$(bat-config --libs)
OBJ = obj/GerdaFitter.o obj/gerda-fitter.o obj/json.gch obj/utils.gch

all : dirs | gerda-fitter $(OBJ)

gerda-fitter : $(OBJ)
	$(CXX) $(FLAGS) -o $@ obj/GerdaFitter.o obj/gerda-fitter.o $(LIBS)

obj/gerda-fitter.o : gerda-fitter.cc json.hpp
	$(CXX) $(FLAGS) -c -o $@ $< $(LIBS)

obj/GerdaFitter.o : GerdaFitter.cc GerdaFitter.hh utils.hpp json.hpp
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
