
PLUGINDIR=$(HOME)/.vmdplugins

ROOT=         cgmap
VERSION=      1.1

EXE=          lib$(ROOT).so
SRC =         cgmap.cpp
OBJ =         $(SRC:.cpp=.o)

CC=           g++ 
CCFLAGS=      -O3 -march=core2 -ffast-math
SHFLAGS=      -fPIC
SHLIBFLAGS=   -shared -m64 

TCL_INC=-I/usr/include/tcl8.5
TCL_STUB_LIB= -ltcl8.5

CGMAP_INC = -DPACKAGE_NAME=\"$(ROOT)\" -DPACKAGE_VERSION=\"$(VERSION)\"

include Makefile.base

all: shlib 

clean:
	rm -f *.o *.so 

install:
	mkdir -p $(PLUGINDIR)/cgtools$(VERSION) && \
		rsync -av ./libcgmap.so $(PLUGINDIR)/cgtools$(VERSION)/.
