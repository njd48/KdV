
#ifndef CARRAY_H
#define CARRAY_H

#include <iostream> 
#include <cassert>
#include <cmath>
#include <complex>
    #define complex_t std::complex<double>
#include <fftw3.h>


class Carray {
private:
        bool isAlloced = false;
        size_t     N;
        complex_t  invN;
        complex_t* z;
        fftw_plan  planf;
        fftw_plan  plani;
public:
        Carray();
        Carray( size_t );
        ~Carray();
        void alloc( size_t );

        complex_t&       operator[]( size_t pos )       { return z[pos]; };
        const complex_t& operator[]( size_t pos ) const { return z[pos]; };

        size_t size(){ return N; };

        void copyTo(         Carray& x );
        void copyFrom( const Carray& x );
        void add(      const Carray& x );
        void add(      const Carray& x, const Carray& y );
        void subtr(    const Carray& x );
        void mult(     const Carray& x );
        void mult(     const Carray& x, const Carray& y );
        void smult(    complex_t a );
        void smult(    double a );
        void addSmult( complex_t a, const Carray&  x);
        void addSmult( double a, const Carray&  x);
        void fft();
        void ifft();

        void disp();
};

Carray::Carray(  ) {}

Carray::Carray( size_t n ) {
    N     = n;
    invN  = (complex_t) ( 1.0 / (double)(N) );
    z     = new complex_t[N];
    planf = fftw_plan_dft_1d( N, reinterpret_cast<fftw_complex*>(z), reinterpret_cast<fftw_complex*>(z), FFTW_FORWARD,  FFTW_ESTIMATE);
    plani = fftw_plan_dft_1d( N, reinterpret_cast<fftw_complex*>(z), reinterpret_cast<fftw_complex*>(z), FFTW_BACKWARD, FFTW_ESTIMATE);
    isAlloced = true;
}

Carray::~Carray(  ) {
    if (isAlloced) {
        delete[] z;
        fftw_destroy_plan(planf);
        fftw_destroy_plan(plani); 
    }
}

void Carray::alloc( size_t n ) {
    N     = n;
    invN  = (complex_t) ( 1.0 / (double)(N) );
    z     = new complex_t[N];
    planf = fftw_plan_dft_1d( N, reinterpret_cast<fftw_complex*>(z), reinterpret_cast<fftw_complex*>(z), FFTW_FORWARD,  FFTW_ESTIMATE);
    plani = fftw_plan_dft_1d( N, reinterpret_cast<fftw_complex*>(z), reinterpret_cast<fftw_complex*>(z), FFTW_BACKWARD, FFTW_ESTIMATE);
    isAlloced = true;
}

void Carray::copyFrom( const Carray& x ){
    assert( N == x.N );
    for ( size_t i = 0; i < N ; i++ ){
        z[i] = x.z[i];
    }
}

void Carray::copyTo( Carray& x ){
    assert( N == x.N );
    for ( size_t i = 0; i < N ; i++ ){
        x.z[i] = z[i];
    }
}


void Carray::add(  const Carray& x ) {
    assert( N == x.N );
    for ( size_t i = 0; i < N ; i++ ){
        z[i] += x.z[i];
    }
}

void Carray::add(  const Carray& x , const Carray& y ) {
    assert( N == x.N );
    assert( N == y.N );
    for ( size_t i = 0; i < N ; i++ ){
        z[i] = x.z[i] + y.z[i];
    }
}

void Carray::subtr( const  Carray& x ) {
    assert( N == x.N );
    for ( size_t i = 0; i < N ; i++ ){
        z[i] -= x.z[i];
    }
}

void Carray::mult( const Carray& x ) {
    assert( N == x.N );
    for ( size_t i = 0; i < N ; i++ ){
        z[i] *= x.z[i];
    }
}

void Carray::mult( const Carray& x, const Carray& y ) {
    assert( N == x.N );
    for ( size_t i = 0; i < N ; i++ ){
        z[i] = x.z[i] * y.z[i];
    }
}

void Carray::smult(  complex_t a ) {
    for ( size_t i = 0; i < N ; i++ ){
        z[i] *= a;
    }
}

void Carray::smult(  double a ) {
    for ( size_t i = 0; i < N ; i++ ){
        z[i] *= a;
    }
}

void Carray::addSmult( complex_t a, const Carray&  x){
    for ( size_t i = 0; i < N ; i++ ){
        z[i] += a*x[i];
    }
}

void Carray::addSmult(    double a, const Carray&  x){
    for ( size_t i = 0; i < N ; i++ ){
        z[i] += a*x[i];
    }
}

void Carray::fft(){
    fftw_execute(planf);
}

void Carray::ifft(){
    fftw_execute(plani);
    smult( this->invN );
}

void Carray::disp(  ) {
    std::cout << "[ ";
    for( size_t n = 0; n<N ; n++ ) {
        std::cout  << z[n] << ' ';
    }
    std::cout << "]\n";
}


#endif
