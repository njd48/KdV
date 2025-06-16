
#ifndef KDV_FSTREAM_H
#define KDV_FSTREAM_H

#include <iostream>
#include <fstream>
#include "Carray.cpp"
#include <complex>
    #define complex_t std::complex<double>

std::ofstream write_file_header( std::string f, int simIdno, size_t N ) {

    std::ofstream fid ( f );

    fid << "KdV simulation \n";
    fid << "Author: Nicholas J. Dubicki \n";
    fid << "\n";
    fid << "Solve the KdV equation on the interval x in [ 0, 2*pi )\n";
    fid << "simIDno: " << simIdno << '\n';
    fid << "N:       " << N << '\n';
    fid << "first column is time t\n";
    fid << "--------------------------------------------\n";


    return fid;

}

void write_data_line( std::ofstream& fid, size_t N, double t, const Carray& U ) {

    fid << t ;
    for( size_t i = 0 ; i<N ; i++ ){
        fid << ',' << U[i].real();
    }
    fid << '\n';
}


#endif