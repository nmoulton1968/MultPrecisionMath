CXXFLAGS =	-O3 -g -Wall

pi:	pi.o mpim.o
	$(CXX) -o pi.exe pi.o mpim.o

e:	e.o mpim.o
	$(CXX) -o e.exe e.o mpim.o

pi.o :	pi.cpp
	$(CXX) -c pi.cpp $(CXXFLAGS)

e.o :	e.cpp
	$(CXX) -c e.cpp $(CXXFLAGS)

mpim.o :	mpim.cpp mpim.h
	$(CXX) -c mpim.cpp $(CXXFLAGS)

clean:
	rm -f -v *.o *.orig pi.exe e.exe
