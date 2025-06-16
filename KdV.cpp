

#ifndef KDV_H
#define KDV_H

#include <iostream>
#include <fstream>
#include <cstddef>
#include <ctime>
#include <cmath>
#include <complex>
#include "fftw3.h"
    #define complex_t std::complex<double>

#include "Carray.cpp"
#include "KdV_fncs.cpp"
#include "KdV_fstream.cpp"


class KdV_env {

private:

    // --- Env size ---------------
    size_t N;  // = DOMAINSIZE;
    
    double xL;
    double xR;

    // --- Principle Simulation Variables ---  
    double t = 0.0; 

    Carray Udata0;  // Data addresses
    Carray Udata1;
    Carray Udata2;
    Carray temp1;
    Carray ikvals;

    Carray* U0; // Data Aliases
    Carray* U1;
    Carray* U2;

    // --- fftw plans -----------------
    //fftw_plan  p0, p1, p2, p0i, p1i, p2i;
    //fftw_plan* P0i, P1i, P2i;

    void swap_alias();

    // --- State vars -----------------
    bool initialized       = false; // initial condition has been applied
    bool inFrequencyDomain = true;  // true if main computations occur in fft space
                                    // necessary for file writes

    // --- Derived parameters ---------
    double dx;

    // --- File Properties ------------
    std::string outDir;
    std::string fname;
    int simIDno;
    std::ofstream fid;
    bool fHasHeader = false;
    bool fisOpen    = false;

    void openfile(  ){ 
        std::string f = outDir + '/' + fname + ".txt";
        if ( fHasHeader ) {
            std::ofstream fid(f);
        } else {
            fid =  write_file_header( f, simIDno, N ); 
        }
        fisOpen = true;
    };
    void closefile( ){ 
        fid.close(); 
        fisOpen = false; 
    }

public:
    // --- Public parameters ------
    double dt      = 0.0;
    double t_final = 0.0;
    double t_write = 0.0;

    // --- constr destr -------
    KdV_env( size_t );
    ~KdV_env();

    // --- Retrievals ---------
    double time() { return t; }

    // --- Ev -----------------
    void set_out_dir( std::string );
    void set_filename( std::string );

    void set_initial();
    void write2file();
    void run();

    void dispDiag();
    void dispDiagifft();

};

//------------------------------------------------------------------------------
// constr and destr
//------------------------------------------------------------------------------

KdV_env::KdV_env( size_t n ){
    
    N = n;

    dx  = (2.*M_PI) / (double)(N); // Set domain variables
    xL  = 0;
    xR  = 2.*M_PI - dx;

    Udata0.alloc(N); // Simulation Data
    Udata1.alloc(N);
    Udata2.alloc(N);
    temp1.alloc(N);
    ikvals.alloc(N);
    KdVcalcs::setiKvals( N, ikvals );

    U0 = &Udata0; // Data Aliases
    U1 = &Udata1;
    U2 = &Udata2;

    outDir = "simData";

}

KdV_env::~KdV_env(){

    if (fisOpen) {
        fid.close();
        fisOpen = false;
    }
    
}

//------------------------------------------------------------------------------
// private methods
//------------------------------------------------------------------------------

void KdV_env::swap_alias() {
    Carray* uuu = U0;
    U0 = U1;
    U1 = U2;    
    U2 = uuu;
}

//------------------------------------------------------------------------------
// public methods
//------------------------------------------------------------------------------

void KdV_env::set_out_dir( std::string s ) {
    outDir   = s;
}

void KdV_env::set_filename( std::string s ) {
    fname   = s;
    simIDno = rand()%10000 * ((int)std::clock()); 
}


void KdV_env::set_initial() {

    using namespace std::complex_literals;
    complex_t x = 0.0 + 0.0*1i;

    for ( size_t n = 0; n<N ; n++ ){
        Udata2[n]= exp( -pow( 5.0*(x-M_PI), 2 )  );
        x += dx;
    }

    dispDiag();

    Udata2.fft();

    dispDiag();
    dispDiagifft();

    // Numerical initial condition for 1st timestep
    // go backward with euler step
    Udata1.copyFrom( Udata2 );
    Carray* u = &Udata0;
    Carray* v = &Udata1;
    Carray* w;
    for ( int tt = 0 ; tt < 4 ; tt++ ) {
        KdVcalcs::reverseEulerTimestep( dt/4.0, dx, ikvals, *u, *v, temp1 );   
        w = u;
        u = v;
        v = w;
    }
    // copy( N, *u, Udata1 ); should not be needed if above loop has even number of iterations  

    initialized = true;
}

void KdV_env::write2file() {
    // std::cout<< "here writing at time, t = " << t << '\n';

    if (inFrequencyDomain) {
        temp1.copyFrom( *U2 );
        temp1.ifft();
        write_data_line( fid, N, t, temp1 );
    } else {
        write_data_line( fid, N, t, *U2 );
    }
    

}

void KdV_env::run() {

    bool eflag = false;
    bool write = true;
    if ( !initialized ) {
        std::cout << "error: system initial condition was not set. cannot run.\n";
        eflag = true;
    } 
    if ( dt <= 0 ) {
        std::cout << "error: invalid timestep. cannot run.\n";
        eflag = true;
    } 
    if ( t_final <= 0 ) {
        std::cout << "error: specify finite t_final. cannot run.\n";
        eflag = true;
    }
    if ( t_write <=0 ) {
        std::cout << "warning: invalid write period, t_write.  output will not be saved. \n";
        write = false;
    }
    if ( t_write < dt ) {
        std::cout << "warning: t_write < dt, replaced with t_write = dt. \n";
        t_write = dt;
    }

    if (eflag) { 
        return; 
    }
    //----------------------------------------------------------------------

    std::cout<< "computing...\n";

    size_t P   = (size_t)( t_final / dt );

    if (write) {

        if ( fname.length() < 1 ) { fname = "u"; }

        openfile();

        size_t P_w = (size_t)( t_write / dt );

        size_t n_writes = P/P_w;

        write2file();

        //int timeIndex = 0;

        for ( size_t n =0 ; n < n_writes ; n++ ) {
            for ( size_t p = 0; p < P_w ; p++ ){

                swap_alias();
                KdVcalcs::timestep( dt, dx, ikvals, *U0, *U1, *U2, temp1 );                
                t += dt;
                //timeIndex++;
                //std::cout << timeIndex << '\n';
            }
            write2file();
        }

        if ( n_writes*P_w < P ) {
           // std::cout << "amhere?\n";
            for ( size_t p = 0; p < P-n_writes*P_w  ; p++ ){
                swap_alias();
                KdVcalcs::timestep( dt, dx, ikvals, *U0, *U1, *U2, temp1 );                
                t += dt;
                //timeIndex++;
              //  std::cout << timeIndex << '\n';
            }
            write2file();
        }
       // std::cout << "expected final time index: " << P << '\n';

    } else {
        for ( size_t p = 0; p < P ; p++ ) 
        {
            swap_alias();
            KdVcalcs::timestep( dt, dx, ikvals, *U0, *U1, *U2, temp1 );                
            t += dt;
        }
    }

}

void KdV_env::dispDiag(){
    U2->disp();
}
void KdV_env::dispDiagifft(){
    U2->copyTo(temp1);
    temp1.ifft();
    temp1.disp();
}


#endif
