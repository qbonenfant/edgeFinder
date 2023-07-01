// SeqAn includes
#include <seqan/arg_parse.h>

// STD C++ includes
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <unordered_set>
#include <fstream>

// custom includes
#include "params.h"
#include "types.h"
#include "approx_ressearch.h"
#include "utils.h"


//#define EF_DEBUG 1 
#if EF_DEBUG and not NDEBUG
#define DEBUG(x) std::cout << x << "\n"
#else
#define DEBUG(x)
#endif

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
//-------------------------  PRINT FUNCTIONS --------------------------------//

// Shortcut to print things with timestamp

// TODO: FIND OUT WHY THE LINKER REFUSE TO FETCH THAT SPECIFIC FUNCTION FROM
// UTIL.H BUT WILL WORKS FINE WITH EVEYTHING ELSE....
// EDIT: Apprently it is because of the template... Well i'll find a way
const auto boot_time = std::chrono::steady_clock::now();
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}

// Simply format and print both elements from each pair of a pair vector.
template<typename TPrintType>
void print_pair_vector(TPrintType my_map){

    for(auto element: my_map){
        std::cout << "(" << element.first << "," << element.second << "), ";
    }
    std::cout << "\n" ;
}

// Similar but with result maps (read id + list of pairs of positions)
template<typename TPrintType>
void print_result_map(TPrintType my_map){
    for(auto element: my_map){
        std::cout << element.first << "\n";
        print_pair_vector(element.second);
    }
    std::cout << "\n" ;
}
//---------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// DEPRECATED This function need to be elsewhere.
// void split_edge_on_degree(std::string edge_file, node_degree_t node_degree){
//     print("SPLITTIN EDGE FILE USING NODE DEGREE");
//     std::ifstream edge_file_stream(edge_file);
//     std::ofstream repr_file(edge_file + "_repr");
//     std::ofstream residual_file(edge_file + "_residual");
//     std::vector<std::string> line_fields;

//     unsigned representant_lines = 0;
//     unsigned residual_lines = 0;
//     if( edge_file_stream.is_open() and repr_file.is_open() and residual_file.is_open()){
//         for( std::string line; getline( edge_file_stream, line );){

//             line_fields = split(line);

//             // Getting the name of current target, sources are always repr
//             if(line_fields.size() > 2){
//                 std::string target = line_fields[2];

//                 // checking if node has degree of 1
//                 if( node_degree[target] == 1){
//                     residual_file << line << "\n";
//                     residual_lines+=1;
//                 }
//                 // If not, export to repr
//                 else{
//                     repr_file << line << "\n";
//                     representant_lines +=1;
//                 }
//             }
//         }
//         print("Exported " + std::to_string(residual_lines) + " residual edges and " +
//               std::to_string(representant_lines) + " representant edges");
//         edge_file_stream.close();
//         repr_file.close();
//         residual_file.close();

//     }
//     else{
//         print("/!\\ \tCOULD NOT SPLIT EDGE FILE: FILES COULD NOT BE OPENED PROPERLY");
//     }

// }






// Search links between reads, which correspond to edges in our graphs.
void find_edges(const seq_id_set_t & ids, const seq_set_t & sequences, index_t & index, std::ofstream & output_file, ef_params & params){
    // Initialising containers
    read2pos_map_t  results;

    // keeping track of already processed reads for "sampling" method
    bool_vector_t processed_reads(length(sequences), false);
    unsigned nb_sequences = length(sequences);

    // keeping trace of number of seeds.
    unsigned total_seed_number = 0;
    unsigned filtered_seed_number = 0;

    // revcomp management
    int is_revcomp = params.rc ? 2 : 1;

    // Processing all sequences
    for(unsigned current_read = 0; current_read < nb_sequences; ++current_read){

        // Print progression
        if(params.v >=2 ){
            if(nb_sequences >=100 and current_read != 0 and  current_read % (nb_sequences/100)  == 0 ){
                print( std::to_string( (current_read*100)/(nb_sequences) ) + "%" );
            }
        }

        // Should we process this read ?
        if(not params.sampling or not processed_reads[current_read]){
            
            // We are processing this read, so do not do it again
            processed_reads[current_read] = true;
            
            for(int revcomp = 0; revcomp < is_revcomp; revcomp ++){
                
                // Clearing containers
                results.clear();
                
                // Ressearch
                find_kmers(index, sequences[current_read], results, current_read, revcomp == 1, params );

                // Filtering edges (LIS, removing mapping below minimim number of common k-mers)
                filter_edges(results, sequences, current_read, revcomp == 1 , params);
                
                // Computing chimera score
                float chimera_score = is_chimera(results, length(sequences[current_read]));
                if(chimera_score < params.ct){
                    //Export to output file
                    export_edges(results, current_read, ids, sequences, revcomp == 1, chimera_score, output_file);
                    
                    // Considering exported reads as processed
                    for(auto el : results){
                        processed_reads[el.first] = true;
                    }
                }
            }
        }
    }

    if(params.v>= 2){
        print("Total number of seeds: " + std::to_string(total_seed_number) );
        print("Filtered number of seeds: " + std::to_string(filtered_seed_number) );
    }

}







///////////////////////////////////////////////////////////////////////////////
//----------  MAIN  ---------------------------------------------------------//
//                                                                           //
int main(int argc, char const ** argv)
{

    //--------------  ARGS PARSING  -------------------------------------------
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("approxCount");

    // 1st position argument: input file
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
        "i", "index", "Path and prefix of the index files",
        seqan::ArgParseArgument::STRING, "index file"));

    addOption(parser, seqan::ArgParseOption(
        "si", "save_index", "Path to the folder to store index files, if wanted.",
        seqan::ArgParseArgument::STRING, "index folder"));

    addOption(parser, seqan::ArgParseOption(
        "o", "outFile", "Path to the output file",
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
        "ks", "kmer-skip", "Skip 'ks' k-mers when ressearching (default: 1 over 3)",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "lc", "low_complexity", "Threshold for low complexity filter",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, seqan::ArgParseOption(
        "kp", "km_prop", "Minimum proportion of common k-mers between two reads.\n\
        Threshold value is computed using the smallest read length, and can not be below nk",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, seqan::ArgParseOption(
        "ct", "chimera_threshold", "Threshold for chimera detection.\n\
        Reads with values above this level will be considered potential chimera.",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));


    addOption(parser, seqan::ArgParseOption(
        "rc", "rev_comp", "Ressearsh each k-mer's rev-comp too, in the case read orientation is mixed. Set to 0 to deactivate.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "lis", "lis_mode", "Choose between spaced LIS (1) or normal LIS(0)",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "s", "star", "Use star method. This is the default implementation.\n\
        Already associated reads will not start new k-mer ressearch. Set to 0 to deactivate",
        seqan::ArgParseArgument::INTEGER, "INT"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;


    isSet(parser, "rev_comp");

    // Creating parameter structure
    ef_params params;

    // Fetching data from parser
    params.set_param_from_args(parser);

    //--------------------- END OF ARG PARSING --------------------------------


    // Creating output file stream
    std::ofstream output_file;
    output_file.open(params.output);

    // Exporting run params to outfile as trace
    params.export_params(output_file);


    // Fasta parsing
    fasta_pair_t fasta;
    
    if(params.v>=1)
        print("PARSING FILE");
    
    parse_fasta(fasta, params);

    if(params.v >= 2){
        print("Number of read parsed:");
        print(length(fasta.first));
    }

    if(params.v>=1)
        print("DONE");


    // Index creation
    if(params.v>=1)
        print("INITIALISING INDEX");
    index_t index(fasta.second);
    create_index(index, params);

    // Ressearch edges and export the results
    find_edges(fasta.first, fasta.second, index, output_file, params);

    // Closing output file
    output_file.close();

    
    if(params.v>=1){
        print("PROGRAM END");
    }

    return 0;
}
//-------------- ##### END ##### --------------