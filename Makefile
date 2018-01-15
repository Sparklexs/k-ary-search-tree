CXX = g++

ifeq ($(DEBUG),1)
    CXXFLAGS = -march=native -mavx2 -std=c++11 -Weffc++ -pedantic -D_GLIBCXX_DEBUG -DDEBUG=1 -ggdb -Wall -Wextra -Wcast-align -Wconversion -Winline
else
    CXXFLAGS = -march=native -mavx2 -std=c++11 -Weffc++ -DNDEBUG=1 -pedantic -O3 -Wall -Wextra -Winline -Wcast-align -Wconversion
endif

HEADERS = $(shell ls include/*.h*)

OBJS =	 intersection.o# k-ary_search_tree.o k-ary_search_tree_trial.o

LIBS =

TARGET =	benchintersection

all: $(TARGET)
	@echo "finish!"

intersection.o: src/intersection.cpp include/common.h
	$(CXX) $(CXXFLAGS) -Iinclude -c src/intersection.cpp
	
#k-ary_search_tree.o: intersection.o $(HEADERS)
#	$(CXX) $(CXXFLAGS) -Iinclude -o k-ary_search_tree.o -c include/k-ary_search_tree.hpp
#	
#k-ary_search_tree_trial.o: intersection.o $(HEADERS)
#	$(CXX) $(CXXFLAGS) -Iinclude -o k-ary_search_tree_trial.o -c include/k-ary_search_tree_trial.hpp

benchintersection: $(HEADERS) $(OBJS) src/benchintersection.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -o benchintersection src/benchintersection.cpp intersection.o
	
#$(TARGET):	$(OBJS) $(HEADERS)
#	$(CXX) $(CXXFLAGS) -Iinclude -o $(TARGET) src/benchintersection.cpp $(OBJS) $(LIBS)


clean:
	rm -f $(OBJS) $(TARGET)
