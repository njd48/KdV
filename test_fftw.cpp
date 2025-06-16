
#include <iostream>
#include <cmath>
#include <complex>
#include <fftw3.h>

#include "Carray.cpp"

void mult( int N, complex_t* x, complex_t* y, complex_t* z) {
    for( int i = 0; i < N ; i++ ) {
        z[i] = x[i]*y[i];
    }
}

void smult( int N, complex_t a, complex_t* x, complex_t* z) {
    for( int i = 0; i < N ; i++ ) {
        z[i] = a*x[i];
    }
}

void dispV( int N , complex_t* U ) {
    std::cout << "[ ";
    for( int n = 0; n<N ; n++ ) {
        std::cout  << U[n] << ' ';
    }
    std::cout << "]\n";
}

void setVals( int N , complex_t* U ) {
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

void setKvals( int N, complex_t* kvals ) {

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

    complex_t* k = new complex_t[N];
    complex_t* u = new complex_t[N];
    complex_t* v = new complex_t[N];

    fftw_plan p    = fftw_plan_dft_1d( N, reinterpret_cast<fftw_complex*>(u), reinterpret_cast<fftw_complex*>(v), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan pinv = fftw_plan_dft_1d( N, reinterpret_cast<fftw_complex*>(v), reinterpret_cast<fftw_complex*>(u), FFTW_BACKWARD, FFTW_ESTIMATE);

    std::cout << '\n';
    std::cout << "Test 1: Try to use fftw_plans to transform u into v and reverse...\n";

    setVals(  N, u );
    dispV(    N, u );
    fftw_execute( p );
    dispV(    N, v );
    fftw_execute( pinv );
    smult( N, 1/(double)(N) , u, u );
    dispV(    N, u );

    std::cout << '\n';
    std::cout << "Try to FT in-place, and alloc k vals and take derivative to arrive at cosine...\n";

    p    = fftw_plan_dft_1d( N, reinterpret_cast<fftw_complex*>(u), reinterpret_cast<fftw_complex*>(u), FFTW_FORWARD, FFTW_ESTIMATE);
    pinv = fftw_plan_dft_1d( N, reinterpret_cast<fftw_complex*>(u), reinterpret_cast<fftw_complex*>(u), FFTW_BACKWARD, FFTW_ESTIMATE);

   
    setKvals( N, k);
    using namespace std::complex_literals;
    smult( N, 1.0*1i, k , k  );
    setVals(  N, u );
    std::cout << "u    = "; dispV(    N, u );
    fftw_execute( p );
    std::cout << "uk   = ";  dispV(    N, u );
    std::cout << "k    = ";  dispV(    N, k );
    mult( N, k, u, u);
    std::cout << "k*uk = ";  dispV(    N, u );
    fftw_execute( pinv );
    smult( N, 1/(double)(N) , u, u );
    std::cout << "v    = "; dispV(    N, u );
    /**/

    fftw_destroy_plan(p); 
    fftw_destroy_plan(pinv); 

    std::cout << "success.\n";

    delete[] k;
    delete[] u;
    delete[] v;
}