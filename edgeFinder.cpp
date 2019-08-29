// SeqAn includes
#include <seqan/arg_parse.h>

// STD C++ includes
#include <iostream>
#include <algorithm> 
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <unordered_set>

// custom includes
#include "approx_ressearch.h"
#include "utils.h"

using namespace seqan;

void export_edges(read2pos_map_t & results, read_id_t current_read_id, const seq_id_set_t & realIds, const seq_set_t & sequences, pos_vector_t & processed_reads, isoform_map_t & iso_map, uint8_t nk,  bool reverse,  std::ofstream & output_file){
    
    for( auto it=results.begin(); it != results.end(); ++it){
        
        // Avoid printing the same id, and results that are too short
        if(it->first != current_read_id and it->second.size() >= nk ){
            processed_reads.emplace_back(it->first);
            output_file << realIds[current_read_id];
            output_file << "\t" << length(sequences[current_read_id]);
            output_file << "\t" << realIds[it->first];
            output_file << "\t" << length(sequences[it->first]);
            output_file << "\t" << reverse;
            output_file << "\t" << iso_map[it->first];

            // Exporting seeds.
            for(auto pos: it->second){
                    output_file << "\t" << pos.first << "," << pos.second;    
            }
            output_file << "\n";
        }
    }
}

isoform_map_t filter_edges(read2pos_map_t & results, const seq_set_t & sequences, unsigned current_read, uint8_t k, uint8_t ks, bool rev_comp, float max_diff_rate, float min_cover){
    isoform_map_t iso_map;
    for(auto read_results: results){
        // If we work with antiparallel strands
        if(rev_comp){
            std::reverse( read_results.second.begin(), read_results.second.end());
        }
        // Keeping longest colinear streak of seeds
        read_results.second = LIS_Pair(read_results.second);
        
        // Checking if both reads are the same isoforms
        // TODO: Refactor is_iso to take the two sequences directly, and derive a "proper coverage" from  k-mer that have been explored
        unsigned l1 = length(sequences[current_read]);
        unsigned l2 = length(sequences[read_results.first]);
        iso_map[read_results.first] = is_iso(read_results.second, l1, l2, k, ks, max_diff_rate, min_cover);


    }
    return(iso_map);
}


// Search links between reads, which correspond to edges in our graphs.
void find_edges(const seq_id_set_t & ids, const seq_set_t & sequences, index_t & index, std::ofstream & output_file, uint8_t k, uint8_t ks, uint8_t nk, uint8_t nb_thread, double lc, bool rc, bool sampling, double mc, double mdr, uint8_t v){

    // Initialising containers
    read2pos_map_t  results;
    isoform_map_t iso_map;
    // keeping track of already processed reads for "sampling" method
    pos_vector_t processed_reads;
    processed_reads.reserve(length(sequences));

    // Processing all sequences
    for(unsigned current_read = 0; current_read < length(sequences); ++current_read){
        
        if(not(sampling) or std::find(processed_reads.begin(), processed_reads.end(), current_read) == processed_reads.end()){

            // Clearing containers
            results.clear();
            iso_map.clear();

            // Forward ressearch
            find_kmers(index, sequences[current_read], results, current_read, k, ks, lc, false);
            
            // Filtering seeds and keeping track of same isoforms
            iso_map = filter_edges(results, sequences, current_read,  k, ks , false, mdr, mc);
            
            //Export to output file
            export_edges(results, current_read, ids, sequences, processed_reads, iso_map, nk, false, output_file);

            if(rc){
                // Clearing containers
                results.clear();
                iso_map.clear();

                // Reverse ressearch
                find_kmers(index, sequences[current_read], results, current_read, k, ks, lc, true);

                // Filtering seeds and keeping track of same isoforms
                iso_map = filter_edges(results, sequences, current_read,  k, ks , true, mdr, mc);
                
                //Export to output file
                export_edges(results, current_read, ids, sequences, processed_reads, iso_map, nk, true, output_file);
            }
        }
    }
}



//--------------  ARGS PARSING  --------------

int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("approxCount");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "input filename"));


    // Describe how the run should go: nb thread, verbosity, output, use an index ?..
    addSection(parser, "Global run parameters");

    addOption(parser, seqan::ArgParseOption(
        "nt", "nb_thread", "Number of thread to work with",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "v", "verbosity", "Level of details printed out (0: nothing, 1: some, 2: a lot)",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "i", "index", "path and prefix of the index files",
        seqan::ArgParseArgument::STRING, "index file"));

    addOption(parser, seqan::ArgParseOption(
        "o", "outFile", "path to the output file",
        seqan::ArgParseArgument::STRING, "output file"));


    // Parameters for the algorithm
    addSection(parser, "Ressearch parameters");    

    addOption(parser, seqan::ArgParseOption(
        "k", "k", "Size of the kmers",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "nk", "nb_kmer", "Minimum number of common k-mer needed to create an edge between two reads.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "ks", "kmer-skip", "Limit ressearch to 1/ks k-mers in the read",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "lc", "low_complexity", "Threshold for low complexity",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, seqan::ArgParseOption(
        "mc", "min_cover", "Minimal seed cover between two reads to consider them as coming from the same isoform",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, seqan::ArgParseOption(
        "mdr", "max_difference_rate", "Maximum offset ratio between two pairs of seeds on reads that come from the same isoform",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, seqan::ArgParseOption(
        "rc", "rev_comp", "Ressearsh each k-mer's rev-comp too, in the case read orientation is mixed. Set to 0 to deactivate.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "s", "sampling", "Use read sampling / no reprocess method. Set to 0 to deactivate",
        seqan::ArgParseArgument::INTEGER, "INT"));


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Setting default values ---------------------------------------------------------
    std::string output = "out.edges";     // output file                               
    std::string index_file = "";  // index prefix, if any                              
    unsigned nb_thread = 4;       // default number of thread (4)                      
    //unsigned nbErr = 2;         // default number of errors EDIT: CAN NOT BE ASSIGNED
    unsigned v = 1;               // verobisty, default = 1                            
    unsigned k = 16;              // kmerSize, default = 16                            
    unsigned ks = 3;              // k-mer skip size, default = 3                      
    unsigned nk= 10;              // Minimum common k-mer, default = 10                 
    double   lc= 1.25;            // "dust2" low complexity threshold, default = 1.25
    double   mc= 0.76;  
    double   mdr= 1.35;           
    bool     rc =  true;          // Rev Comp research ? True by default               
    bool sampling =  true;        // Sampling method   ? True by default                 
    // --------------------------------------------------------------------------------


    // ASSIGNING VALUES FROM PARSER ---------------------------------------------------
    //
    // Boolean values
    // checking if revcomp is activated
    if( isSet(parser, "rev_comp") ){
        int val;
        getOptionValue(val, parser, "rc");
        rc = val!=0;
    }
    
    // checking if sampling is activated
    if( isSet(parser,"sampling") ){
        int val;
        getOptionValue(val, parser, "sampling");
        sampling = val!=0;
    }
    
    // Numerical values
    getOptionValue(nb_thread, parser, "nt");
    getOptionValue(v, parser, "v");
    getOptionValue(k, parser, "k");
    getOptionValue(nk, parser, "nk");
    getOptionValue(ks, parser, "ks");
    getOptionValue(lc, parser, "lc");
    getOptionValue(mc, parser, "mc");
    getOptionValue(mdr, parser, "mdr");
    getOptionValue(index_file, parser, "i");
    
    // Input file (fasta only)
    std::string fasta_file;
    setValidValues(parser, 0, "FASTA fa");
    getArgumentValue(fasta_file, parser, 0);
    
    // Output file 
    getOptionValue(output, parser, "o");
    std::ofstream output_file;
    output_file.open (output);


    // Exporting run parameters to output file
    output_file << "file:" << fasta_file;  
    output_file << " nb_thread:" << nb_thread;
    output_file << " k:" << k;
    output_file << " nk:" << nk;
    output_file << " kmer_skipped:" << ks;
    output_file << " lc:" << lc;
    // adjusting threshold to current kmer size
    lc = adjust_threshold( lc, 16, k );
    // display new threshold
    output_file << " adusted_lc:" << lc;
    output_file << " rev_comp:" << rc;
    output_file << " sampling:" << sampling;
    output_file << " index:" << index_file;
    output_file <<std::endl;
    output_file <<  "output:" << output << std::endl;

    // PROGRAM STARTING POINT

    // Fasta parsing
    fasta_pair fasta;
    parse_fasta(fasta_file, fasta, v);

    // Index creation
    
    if(v>=1)
        print("INITIALISING INDEX");
    index_t index(fasta.second);
    create_index(index_file, index, v);

    // Ressearch edges and export the results
    find_edges(fasta.first, fasta.second, index, output_file, k, ks, nk, nb_thread, lc, rc, sampling, mc, mdr, v);
    
    
    // closing
    output_file.close();
    if(v>=1){
        print("PROGRAM END");
    }

    return 0;
}
//-------------- ##### END ##### --------------