CXX:=g++-3.3
MATLAB:=/Applications/_IKW/Matlab74/bin/
MEX:=$(MATLAB)mex -v
LDFLAGS:=-L$(MATLAB)mac -lm libedfapi.a -lstdc++
SOURCE:=../../src/
INCLUDE:=../../include/
default: edfread.mexmac

all: edfread.mexmac static

edfread.mexmac:
	$(MEX) -f mexopts.sh $(LDFLAGS) -DMEX -I$(SOURCE) -I$(INCLUDE) $(SOURCE)edfread.cpp

static:
	$(CXX) -g $(SOURCE)edfread.cpp -o edfread.mac.static -I$(SOURCE) -I$(INCLUDE) $(LDFLAGS) -lgcc


clean:
	rm -f edfread.mexmac
	rm -f edfread.mac.static
