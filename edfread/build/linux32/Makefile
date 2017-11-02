CXX:=g++-3.3
#MATLAB:=/home/jss/bin/ml/bin/
MATLAB:=/opt/matlab/bin/
#MATLAB:=/Applications/_IKW/Matlab704/bin/
MEX:=$(MATLAB)mex -v
LDFLAGS:=-L$(MATLAB)glnx86 -lm
SOURCE:=../../src/
INCLUDE:=../../include/
default: edfread.mexglx

all: edfread.mexglx static

edfread.mexglx:
	$(MEX) -f mexopts.sh $(LDFLAGS) -DMEX -I$(INCLUDE) -I$(SOURCE) $(SOURCE)edfread.cpp libedfapi.a

static:
	$(CXX) -g $(SOURCE)edfread.cpp libedfapi.a -o edfread.x86.static -I$(INCLUDE) -I$(SOURCE) $(LDFLAGS) -lgcc 

clean:
	rm -f edfread.mexglx
	rm -f edfread.x86.static
