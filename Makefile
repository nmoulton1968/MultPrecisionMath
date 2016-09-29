CXXFLAGS =	-O3 -g -Wall
OBJS = 		Pi.o mpim.o
LIBS =

TARGET = 	Pi.exe

$(TARGET):	$(OBJS) mpim.h
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

mpim.o :	mpim.cpp mpim.h
	$(CXX) -c mpim.cpp $(CXXFLAGS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
