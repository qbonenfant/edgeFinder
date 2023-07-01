#pragma once

#include "version.h"
#include "types.h"
#include <fstream>
#include <seqan/arg_parse.h>

// Adjust the low complexity threshold to k-mer size, usin a k_old as base
float adjust_threshold(float c_old, uint8_t k_old, uint8_t k_new );

//-------------- ARGS STORING  --------------
struct ef_params{

    // Setting default values
    std::string fasta_file;               // input file, can be fasta DNA
    std::string output = "out.edges";     // output file
    std::string index_file   = "";        // index prefix, if any
    std::string index_folder = "";        // saving the index here if specified
    unsigned v  = 1;                      // verobisty, default = 1
    unsigned k  = 15;                     // kmerSize, default = 15
    unsigned ks = 3;                      // k-mer skip size, default = 3
    unsigned nk = 10;                     // Minimum common k-mer, default = 10
    double   lc = 1.00;                   // "dust2" low complexity threshold, default = 1.00
    bool     rc =  true;                  // Rev Comp research ? True by default
    float    kp = 0.0;                    // Minimum k-mer prop. If 0, only nk is used
    float    ct = 0.75;                   // Chimera score threshold
    unsigned nb_thread  =  4;             // default number of thread (4) // INACTIVE FOR NOW
    bool     lis_mode =  true;            // LIS mode, Spaced (1) by default
    bool     sampling   =  true;          // Star method   ? True by default

    void export_params(std::ofstream & output_file);
    void set_param_from_args( seqan::ArgumentParser & parser);


};