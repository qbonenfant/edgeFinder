#include "params.h"


// adjust threshold value to kmer size
float adjust_threshold(float c_old, uint8_t k_old, uint8_t k_new ){
    float c_new = c_old * float( std::pow(k_new - 2 + 1,2) /  std::pow(k_old - 2 + 1,2));
    return(c_new);
}

///////////////////////////////////////////////////////////////////////////////
//----------  PARAMETER STRUCT METHOD DECLARATION  --------------------------//
// This may not be the proper way to do this...                              //
// Struct declaration in type.h

// Export args to outfile
void ef_params::export_params(std::ofstream & output_file){

    output_file << "file: " << fasta_file;
    output_file << " nb_thread: " << nb_thread;
    output_file << " k: " << k;
    output_file << " nk: " << nk;
    output_file << " kmer_skipped: " << ks;
    output_file << " adusted_lc: " << lc;
    output_file << " rev_comp: " << rc;
    output_file << " sampling: " << sampling;
    output_file << " LIS_mode: " << lis_mode;
    output_file <<std::endl;
    output_file <<  "Program version: " << VERSION;
    output_file << " kmer_min_prop: " << kp;
    output_file << " chimera_threshold: " << ct;
    output_file << " index: " << index_file  <<std::endl;
    output_file <<  "output: " << output << std::endl;
}


//---------ASSIGNING VALUES FROM PARSER ---------
void ef_params::set_param_from_args(seqan::ArgumentParser & parser){    
    // Boolean values
    // checking if revcomp is activated
    if( isSet(parser, "rev_comp") ){
        int val;
        getOptionValue(val, parser, "rc");
        rc = val!=0;
    }

    // checking LIS mode
    if( isSet(parser, "lis_mode") ){
        int val;
        getOptionValue(val, parser, "lis_mode");
        lis_mode = val!=0;
    }

    // checking if sampling is activated
    if( seqan::isSet(parser,"star") ){
        int val;
        getOptionValue(val, parser, "star");
        sampling = val!=0;
    }

    // Numerical values
    getOptionValue(nb_thread, parser, "nt");
    getOptionValue(v, parser, "v");
    getOptionValue(k, parser, "k");
    getOptionValue(nk, parser, "nk");
    getOptionValue(ks, parser, "ks");
    getOptionValue(lc, parser, "lc");
    getOptionValue(kp, parser, "kp");
    getOptionValue(ct, parser, "ct");


    // adjusting threshold to current kmer size
    // based on k = 16 as "standard" 
    lc = adjust_threshold( lc, 16, k );

    // Index file, if any
    getOptionValue(index_file, parser, "i");
    // Index storing folder if you want to keep it.
    getOptionValue(index_folder, parser, "si");

    // Input file
    setValidValues(parser, 0, "FASTA fa");
    getArgumentValue(fasta_file, parser, 0);

    // Output file
    getOptionValue(output, parser, "o");

}
//                                                                           //
//---------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////

