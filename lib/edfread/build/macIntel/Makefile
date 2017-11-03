CXX:=g++
MATLAB:=/Applications/_IKW/Matlab74/bin/
MEX:=$(MATLAB)mex -v
LDFLAGS:=-L$(MATLAB)maci -lm
SOURCE:=../../src/
INCLUDE:=../../include/

default: edfread.mexmaci

all: edfread.mexmaci static

edfread.mexmaci:
	$(MEX) -f mexopts.sh $(LDFLAGS) -DMEX -I$(SOURCE) -I$(INCLUDE) $(SOURCE)edfread.cpp edfapi
	install_name_tool edfread.mexmaci -change  " @executable_path/../Frameworks/edfapi.framework/Versions/A/edfapi" "edfapi"


static:
	$(CXX) -g $(SOURCE)edfread.cpp edfapi -o edfread.maci.static -I$(INCLUDE) -I$(SOURCE) $(LDFLAGS) -lgcc
	install_name_tool edfread.maci.static -change  " @executable_path/../Frameworks/edfapi.framework/Versions/A/edfapi" "edfapi"


clean:
	rm -f edfread.mexmaci
	rm -f edfread.maci.static
