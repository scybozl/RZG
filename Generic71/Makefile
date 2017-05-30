CPPFLAGS=`herwig-config --cppflags`
LDFLAGS=`herwig-config --ldflags`
LDLIBS=`herwig-config --ldlibs`

all: ETScale.h ETScale.cc
	g++ -std=c++11 --shared -fPIC -I . $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) -o ETScale.so ETScale.cc

