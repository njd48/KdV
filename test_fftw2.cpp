
#include <iostream>
#include <cmath>
#include <complex>
#include <fftw3.h>

#include "Carray.cpp"

/*
void dispV( int N , complex_t* U ) {
    std::cout << "[ ";
    for( int n = 0; n<N ; n++ ) {
        std::cout  << U[n] << ' ';
    }
    std::cout << "]\n";
}
*/
void setVals( int N , Carray& U ) {
    using namespace std::complex_literals;
    std::complex<double> x  = 0.0 + 0.0*1i;
    std::complex<double> dx = (2.0*M_PI)/(double)(N) + 0.*1i;
    std::complex<double> y;

    for ( int n = 0; n<N ; n++ ){
        y = sin( x );
        U[n] = y;
        x = x +  dx;
    }
}


void setKvals( int N, Carray& kvals ) {

    for ( int k = 0 ; k < (N+1)/2 ; k++ ) {
        kvals[k]   =  k;
        //kvals[N-k] = -k;
    }
    for ( int k = (N+1)/2 ; k < N ; k++  ) {
        kvals[k] = k-N;
    }

}

/*
     fftw_complex in[N], out[N];
     fftw_plan p;
     ...
     p = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE);
     ...
     fftw_one(p, in, out);
     ...
     fftw_destroy_plan(p);  
*/

int main() {

    int N = 8;

    Carray k(N), u(N);

    setKvals( N, k );

    std::cout << '\n';
    std::cout << "Test 3: Try to use fftw embedded in Carray datastructure...\n";

    setVals( N, u);
    u.disp();
    u.fft();
    u.disp();
    u.ifft();
    u.disp();

    std::cout << "success.\n";

}