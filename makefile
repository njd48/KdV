IDIR=./
CXX=g++
CXXFLAGS =-I$(IDIR) -Wall -Wextra -pedantic-errors -O2 -pthread
FFTW_LIB =-lfftw3

.PHONY: clean

test_fftw: test_fftw.cpp test_fftw2.cpp
	$(CXX) test_fftw.cpp -o test_fftw.exe $(FFTW_LIB) $(CXXFLAGS)  
	$(CXX) test_fftw2.cpp -o test_fftw2.exe $(FFTW_LIB) $(CXXFLAGS)  
	./test_fftw.exe
	./test_fftw2.exe

test_KdV: test_KdV.cpp
	$(CXX) test_KdV.cpp -o test_KdV.exe $(FFTW_LIB) $(CXXFLAGS) 
	./test_KdV.exe

clean:
	rm -f test_fftw.exe test_KdV.exe
