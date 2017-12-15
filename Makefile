CXX = g++

ifeq ($(DEBUG),1)
    CXXFLAGS = -march=native -mavx2 -std=c++11 -Weffc++ -pedantic -D_GLIBCXX_DEBUG -DDEBUG=1 -ggdb -Wall -Wextra -Wcast-align -Wconversion -Winline
else
    CXXFLAGS = -march=native -mavx2 -std=c++11 -Weffc++ -DNDEBUG=1 -pedantic -O3 -Wall -Wextra -Winline -Wcast-align -Wconversion
endif

HEADERS = $(shell ls include/*.h*)

OBJS =		k-ary_search_tree.o intersection.o

LIBS =

TARGET =	k-ary_search_tree

all: $(TARGET)
	@echo "finish!"

intersection.o: src/intersection.cpp include/common.h
	$(CXX) $(CXXFLAGS) -Iinclude -c src/intersection.cpp
	
k-ary_search_tree.o: src/k-ary_search_tree.cpp intersection.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -Iinclude -c src/k-ary_search_tree.cpp
	
$(TARGET):	$(OBJS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)


clean:
	rm -f $(OBJS) $(TARGET)
