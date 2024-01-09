.PHONY: all clean test cpptest 

CXX=         g++
CXXFLAGS=    -g -Wall -Wno-unused-function -std=c++17 -O2
LDFLAGS=     -lz 
SRC=         src
TESTS=       tests
GTEST=       $(TESTS)/googletest/googletest


all: prophasm

test: cpptest


cpptest: prophasmtest
	./prophasmtest


prophasm: $(SRC)/main.cpp $(SRC)/$(wildcard *.cpp *.h *.hpp) src/version.h
	./create-version.sh
	$(CXX) $(CXXFLAGS) $(SRC)/main.cpp -o $@ $(LDFLAGS)


prophasmtest: $(TESTS)/unittest.cpp gtest-all.o $(SRC)/$(wildcard *.cpp *.h *.hpp) $(TESTS)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include $(TESTS)/unittest.cpp gtest-all.o -pthread -o $@ $(LDFLAGS)

gtest-all.o: $(GTEST)/src/gtest-all.cc $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

src/version.h: src/version
	./create-version.sh

clean:
	rm -f prophasm
	rm -f prophasmtest
	rm -r -f ./bin
	rm -f gtest-all.o
	rm -f src/version.h
