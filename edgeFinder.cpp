// SeqAn includes
#include <seqan/arg_parse.h>

// STD C++ includes
#include <iostream>
#include <algorithm> 
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <unordered_set>

// custom includes
#include "approx_ressearch.h"
#include "utils.h"

#define VERSION 0.2

#define DEBUG_FLAG 0
#if DEBUG_FLAG and not NDEBUG
#define DEBUG(x) std::cout << x << "\n"
#else
#define DEBUG(x)
#endif

using namespace seqan;

template<typename TPrintType>
void print_pair_vector(TPrintType my_map){

    for(auto element: my_map){
        std::cout << "(" << element.first << "," << element.second << "), ";
    }
    std::cout << "\n" ;
}

template<typename TPrintType>
void print_result_map(TPrintType my_map){   
    for(auto element: my_map){
        std::cout << element.first << "\n";
        print_pair_vector(element.second);
    }
    std::cout << "\n" ;
}


void export_edges(read2pos_map_t & results, read_id_t current_read, const seq_id_set_t & realIds, const seq_set_t & sequences, pos_vector_t & processed_reads, isoform_map_t & iso_map, node_degree_t & node_degree,  uint8_t nk,  bool reverse, node_type_t & node_types, std::ofstream & output_file){
    
    for( auto it=results.begin(); it != results.end(); ++it){
        
        // Avoid printing the same id, and results that are too short
        if(it->first != current_read and it->second.size() >= nk ){
            
            if(iso_map[it->first]>=2){
                node_types[ std::string(toCString(realIds[it->first]))] = true;
            }
            processed_reads.emplace_back(it->first);
            
            
            output_file << realIds[current_read];
            output_file << "\t" << length(sequences[current_read]);
            output_file << "\t" << realIds[it->first];
            output_file << "\t" << length(sequences[it->first]);
            output_file << "\t" << reverse;
            output_file << "\t" << std::to_string(iso_map[it->first]);

            // updating degree
            node_degree[std::string(toCString(realIds[current_read]))] += 1;
            node_degree[std::string(toCString(realIds[it->first]))] += 1;

            // Exporting seeds.
            for(auto pos: it->second){
                    output_file << "\t" << pos.first << "," << pos.second;    
            }
            output_file << "\n";
        }
    }
}

void split_edge(std::string edge_file, node_type_t node_types){

    std::ifstream edge_file_stream(edge_file);
    std::ofstream repr_file(edge_file + "_repr");
    std::ofstream residual_file(edge_file + "_residual");
    std::vector<std::string> line_fields;

    unsigned representant_lines = 0;
    unsigned residual_lines = 0;
    if( edge_file_stream.is_open() and repr_file.is_open() and residual_file.is_open()){
        for( std::string line; getline( edge_file_stream, line );){

            line_fields = split(line);

            // Getting the name of current target, sources are always repr
            if(line_fields.size() > 2){
                std::string target = line_fields[2];
                
                // checking if edge is between isoforms or not
                if( node_types[target]){
                    residual_file << line << "\n";
                    residual_lines+=1;
                }
                // If not, export to repr
                else{
                    repr_file << line << "\n";
                    representant_lines +=1;
                }
            }
        }

        print("Exported " + std::to_string(residual_lines) + " residual edges and " + std::to_string(representant_lines) + " representant edges");
        edge_file_stream.close();
        repr_file.close();
        residual_file.close();

    }
    else{
        print("/!\\ \tCOULD NOT SPLIT EDGE FILE: FILES COULD NOT BE OPENED PROPERLY");
    }

}

void split_edge_on_degree(std::string edge_file, node_degree_t node_degree){
    print("SPLITTIN EDGE FILE USING NODE DEGREE");
    std::ifstream edge_file_stream(edge_file);
    std::ofstream repr_file(edge_file + "_repr");
    std::ofstream residual_file(edge_file + "_residual");
    std::vector<std::string> line_fields;

    unsigned representant_lines = 0;
    unsigned residual_lines = 0;
    if( edge_file_stream.is_open() and repr_file.is_open() and residual_file.is_open()){
        for( std::string line; getline( edge_file_stream, line );){

            line_fields = split(line);

            // Getting the name of current target, sources are always repr
            if(line_fields.size() > 2){
                std::string target = line_fields[2];
                
                // checking if node has degree of 1
                if( node_degree[target] == 1){
                    residual_file << line << "\n";
                    residual_lines+=1;
                }
                // If not, export to repr
                else{
                    repr_file << line << "\n";
                    representant_lines +=1;
                }
            }
        }
        print("Exported " + std::to_string(residual_lines) + " residual edges and " + std::to_string(representant_lines) + " representant edges");
        edge_file_stream.close();
        repr_file.close();
        residual_file.close();

    }
    else{
        print("/!\\ \tCOULD NOT SPLIT EDGE FILE: FILES COULD NOT BE OPENED PROPERLY");
    }

}


void filter_edges(read2pos_map_t & results, const seq_set_t & sequences, isoform_map_t & iso_map, unsigned current_read, uint8_t k, uint8_t ks, bool rev_comp, float max_diff_rate, float min_cover, bool lis_mode){
    for(auto & read_results: results){
        
        // If we work with antiparallel strands
        if(rev_comp){
            std::reverse( read_results.second.begin(), read_results.second.end());
        }
        // Keeping longest colinear streak of seeds
        // DEV-TEST Changing to spaced LIS
        if(lis_mode){
            read_results.second = spaced_LIS_Pair(read_results.second, k);
        }
        // Or normal LIS, depending on option
        else{
            read_results.second = LIS_Pair(read_results.second);
        }

        
        // Checking if both reads are the same isoforms
        // TODO: Refactor is_iso to take the two sequences directly, and derive a "proper coverage" from  k-mer that have been explored
        unsigned l1 = length(sequences[current_read]);
        unsigned l2 = length(sequences[read_results.first]);
        iso_map[read_results.first] = is_iso(read_results.second, l1, l2, k, ks, max_diff_rate, min_cover);

    }
}

// TODO: FIND OUT WHY THE LINKER REFUSE TO FETCH THAT SPECIFIC FUNCTION FROM UTIL.H BUT WILL WORKS FINE WITH EVEYTHING ELSE....
const auto boot_time = std::chrono::steady_clock::now();
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}

unsigned count_seeds(read2pos_map_t results, unsigned current_read ){
    unsigned nb_seeds = 0;
    for( auto it=results.begin(); it != results.end(); ++it){
        
        // Avoid printing the same id, and results that are too short
        if(it->first != current_read){

            nb_seeds += it->second.size();
        }
    }
    return(nb_seeds);
}

// Search links between reads, which correspond to edges in our graphs.
void find_edges(const seq_id_set_t & ids, const seq_set_t & sequences, index_t & index, std::ofstream & output_file, uint8_t k, uint8_t ks, uint8_t nk, uint8_t nb_thread, double lc, bool rc, bool sampling, double mc, double mdr, node_type_t & node_types, node_degree_t & node_degree , bool lis_mode, uint8_t v){
    // Initialising containers
    read2pos_map_t  results;
    isoform_map_t iso_map;

    // keeping track of already processed reads for "sampling" method
    pos_vector_t processed_reads;
    processed_reads.reserve(length(sequences));
    unsigned nb_sequences = length(sequences);

    // keeping trace of number of seeds.
    unsigned total_seed_number = 0;
    unsigned filtered_seed_number = 0;


    // Processing all sequences
    for(unsigned current_read = 0; current_read < nb_sequences; ++current_read){
        
        // Following progress
        if(v <=2 ){
            if(nb_sequences >=100 and current_read != 0 and  current_read % (nb_sequences/100)  == 0 ){
                print( std::to_string( (current_read*100)/(nb_sequences) ) + "%" );
            }
        }

        if(not(sampling) or std::find(processed_reads.begin(), processed_reads.end(), current_read) == processed_reads.end()){

            processed_reads.emplace_back(current_read);

            // Clearing containers
            results.clear();
            iso_map.clear();
            // Forward ressearch
            find_kmers(index, sequences[current_read], results, current_read, k, ks, lc, false);
            
            // Filtering seeds and keeping track of same isoforms
            filter_edges(results, sequences, iso_map, current_read,  k, ks , false, mdr, mc, lis_mode);
            //Export to output file
            export_edges(results, current_read, ids, sequences, processed_reads, iso_map, node_degree, nk, false, node_types, output_file);

            if(rc){
                // Clearing containers
                results.clear();
                iso_map.clear();

                // Reverse ressearch
                find_kmers(index, sequences[current_read], results, current_read, k, ks, lc, true);
                total_seed_number += count_seeds(results, current_read);
                
                // Filtering seeds and keeping track of same isoforms
                filter_edges(results, sequences, iso_map, current_read,  k, ks , true, mdr, mc, lis_mode);
                filtered_seed_number += count_seeds(results, current_read);

                //Export to output file
                export_edges(results, current_read, ids, sequences, processed_reads, iso_map, node_degree, nk, true, node_types, output_file);
            }
        }
    }
    
    if(v>= 2){
        print("Total number of seeds: " + std::to_string(total_seed_number) );
        print("Filtered number of seeds: " + std::to_string(filtered_seed_number) );
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
        "slis", "spaced_lis", "Choose between spaced LIS (1) or normal LIS(0)",
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
    std::string index_file   = "";  // index prefix, if any                              
    std::string index_folder = "";  // saving the index here if specified
    unsigned nb_thread = 4;       // default number of thread (4)                      
    //unsigned nbErr = 2;         // default number of errors EDIT: CAN NOT BE ASSIGNED
    unsigned v = 1;               // verobisty, default = 1                            
    unsigned k = 16;              // kmerSize, default = 16                            
    unsigned ks = 3;              // k-mer skip size, default = 3                      
    unsigned nk= 10;              // Minimum common k-mer, default = 10                 
    double   lc= 1.00;            // "dust2" low complexity threshold, default = 1.00
    double   mc= 0.75;  
    double   mdr= 1.35;           
    bool     rc =  true;          // Rev Comp research ? True by default               
    bool     spaced_lis =  true;  // Spaced LIS mode   ? True by default
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

    // checking LIS mode
    if( isSet(parser, "spaced_lis") ){
        int val;
        getOptionValue(val, parser, "spaced_lis");
        spaced_lis = val!=0;
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

    // Index file, if any
    getOptionValue(index_file, parser, "i");
    // Index storing folder if you want to keep it.
    getOptionValue(index_folder, parser, "si");
    
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
    output_file << " LIS_mode:" << spaced_lis;
    output_file << " index:" << index_file;
    output_file <<std::endl;
    output_file <<  "output:" << output << std::endl;
    output_file <<  "Program version: " << VERSION << std::endl;


    // PROGRAM STARTING POINT

    // Fasta parsing
    fasta_pair fasta;
    parse_fasta(fasta_file, fasta, v);

    // Index creation
    
    if(v>=1)
        print("INITIALISING INDEX");
    index_t index(fasta.second);
    create_index(index_file, index, index_folder, v);

    // Initit the node type map and degree map
    node_type_t node_types;
    node_degree_t node_degree;
    for(auto read_name: fasta.first){
        node_types.insert( std::pair<std::string, bool>( std::string(toCString(read_name)), false) );
        node_degree.insert( std::pair<std::string, unsigned>(std::string(toCString(read_name)), 0) );
    }

    // Ressearch edges and export the results
    find_edges(fasta.first, fasta.second, index, output_file, k, ks, nk, nb_thread, lc, rc, sampling, mc, mdr, node_types, node_degree, spaced_lis, v);
    
    // Closing output file
    output_file.close();

    // counting number of redisudal nodes
    unsigned residual_node_number = 0;
    unsigned representant_node_number = 0;
    
    //for(auto node: node_types){
    for(auto node: node_degree){
        if(node.second==1){
            residual_node_number++ ;
        }
        else{
            representant_node_number++;
        }
    }
    if(v>=1){
        print("Number of residual nodes: " + std::to_string(residual_node_number));
        print("Number of representant nodes: " + std::to_string(representant_node_number));
        print("Total exported nodes: "+ std::to_string(residual_node_number + representant_node_number));
        print("Residual node ratio: " + std::to_string( double(residual_node_number * 100 ) / (residual_node_number + representant_node_number)) + "%" );
        print("---------------------------------------------------");
        print("Total number of nodes: " + std::to_string(length(fasta.first)));
        print("Number of singletons: " + std::to_string(length(fasta.first)-(residual_node_number + representant_node_number)));
    }
    // Splitting graph
    //split_edge(output, node_types);
    split_edge_on_degree(output, node_degree);
    
    
    if(v>=1){
        print("PROGRAM END");
    }

    return 0;
}
//-------------- ##### END ##### --------------