# Makefile
#
# Author: Luigi Pertoldi - pertoldi@pd.infn.it
# Created: Mon 6 May 2019

FLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++11 -g -O3 -I obj \
        $$(bat-config --cflags)
LIBS  = $$(bat-config --libs)
OBJ = obj/HMixFit.o obj/hmixfit.o obj/json.gch obj/utils.gch

all : dirs | hmixfit $(OBJ)

hmixfit : $(OBJ)
	$(CXX) $(FLAGS) -o $@ obj/HMixFit.o obj/hmixfit.o $(LIBS)

obj/hmixfit.o : hmixfit.cc json.hpp
	$(CXX) $(FLAGS) -c -o $@ $< $(LIBS)

obj/HMixFit.o : HMixFit.cc HMixFit.hh utils.hpp json.hpp
	$(CXX) $(FLAGS) -c -o $@ $< $(LIBS)

obj/%.gch : %.hpp
	$(CXX) $(FLAGS) -c -o $@ $^ $(LIBS)

dirs :
	mkdir -p obj

.PHONY : install clean dirs

clean :
	-rm -rf obj hmixfit

install : hmixfit
	install -d $(PREFIX)/bin
	install $< $(PREFIX)/bin
