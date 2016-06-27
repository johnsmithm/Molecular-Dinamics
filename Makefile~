CXX = g++
CXXFLAGS = -Wall -std=c++11 -pedantic
EXTRA = -ftree-vectorizer-verbose=1 
TARGET =waveguide
HXX= READ.h CG.h  Timer.h

OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) $(HXX) 
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf *.o $(TARGET)
