
#ifndef KDV_FUNCS_H
#define KDV_FUNCS_H

#include <complex>
    #define complex_t std::complex<double>

#include "Carray.cpp"
//#include "arraymeth.cpp"


// =============================================================
//  Forward Declarations
// =============================================================
namespace KdVcalcs {

    //fftw_complex* reinterp( complex_t* u );

    void setiKvals( size_t N , Carray& ikvals );
    void mderivFperiodic( size_t N, double dx, const Carray& U,            Carray& D1  );
    void mderivBperiodic( size_t N, double dx, const Carray& U,            Carray& D1  );
    void mderivCperiodic( size_t N, double dx, const Carray& U,            Carray& D1  );
    void mderivKperiodic( const Carray& ikvals, const Carray& U, Carray& D1);
    void RHS( double dx, const Carray& ikvals, const Carray& U, Carray& rhs );
    void timestep( double dt, double dx, const Carray& ikvals, const Carray& u0, const Carray& u1, Carray& u2, Carray& temp );
    void reverseEulerTimestep( double dt,  double dx, const Carray& ikvals, Carray& u0, Carray& u1, Carray& temp );

};

// =============================================================
//  Functions
// =============================================================

    void KdVcalcs::setiKvals( size_t N , Carray& ikvals ){
        using namespace std::complex_literals;
        for ( size_t k = 0 ; k < (N+1)/2 ; k++ ) {
            ikvals[k]   =  (double)(k) * 1i ;
        }
        for ( size_t k = (N+1)/2 ; k < N ; k++  ) {
            ikvals[k] = (double)(k-N) * 1i ;
        }
    }

    void KdVcalcs::mderivFperiodic( size_t N, double dx, const Carray& U, Carray& D1  )
    {

        //D1[0]   = -1.0/(2.0*dx) * ( U[1]-U[N-1] );
        D1[N-1] = -1.0/(dx) * ( U[0]-U[N-1] );

        for( size_t j = 0 ; j < N ; j++ ) {
            D1[j] = -1.0/(dx) * ( U[j+1] - U[j] );
        }
    }

    void KdVcalcs::mderivBperiodic( size_t N, double dx, const Carray& U, Carray& D1  )
    {

        D1[0]   = -1.0/(dx) * ( U[0]-U[N-1] );
        //D1[N-1] = -1.0/(dx) * ( U[0]-U[N-1] );

        for( size_t j = 1 ; j < N-1 ; j++ ) {
            D1[j] = -1.0/(dx) * ( U[j] - U[j-1] );
        }
    }

    void KdVcalcs::mderivCperiodic( size_t N,  double dx, const Carray& U, Carray& D1   )
    {

        D1[0]   = -1.0/(2.0*dx) * ( U[1]-U[N-1] );
        D1[N-1] = -1.0/(2.0*dx) * ( U[0]-U[N-2] );

        for( size_t j = 1 ; j < N-1 ; j++ ) {
            D1[j] = -1.0/(2.0*dx) * ( U[j+1] - U[j-1] );
        }
    }

    void KdVcalcs::mderivKperiodic(  const Carray& ikvals, const Carray& U, Carray& D1   )
    {
        D1.mult(ikvals, U);
        D1.smult(-1.0);
    }

    void KdVcalcs::RHS( double dx, const Carray& ikvals, const Carray& U, Carray& rhs ) 
    {
        //KdVcalcs::mderivBperiodic( rhs.size(), dx, U, rhs );
        //KdVcalcs::mderivCperiodic( N, dx, U, rhs );
        KdVcalcs::mderivKperiodic( ikvals, U, rhs );
    }

    void KdVcalcs::timestep( double dt, double dx, const Carray& ikvals, const Carray& u0, const Carray& u1, Carray& u2, Carray& temp ){

       
        KdVcalcs::RHS(  dx, ikvals, u1,  temp );

        temp.smult( 2.0*dt );
        u2.add( u0, temp);

    }

    void KdVcalcs::reverseEulerTimestep(  double dt, double dx, const Carray& ikvals, Carray& u0, Carray& u1, Carray& temp ){

       
        KdVcalcs::RHS(  dx, ikvals, u1,  temp );

        temp.smult( -dt );

        u0.add( u1, temp );

    }



#endif