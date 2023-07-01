#pragma once

// required types
#include "types.h"
#include "params.h"
#include <seqan/seq_io.h>



// print template, containing a timer from program start
template<typename TPrintType>
void print(TPrintType text);

// Export edge to outfile
void export_edges(read2pos_map_t & results, read_id_t current_read, const seq_id_set_t & realIds, const seq_set_t & sequences, bool reverse, float score, std::ofstream & output_file);

// Parse the input fasta file
void parse_fasta(fasta_pair_t & fasta, ef_params & params);

// Create the FM-Index
void create_index(index_t & index, ef_params & params);

// Split strings on delimiter
std::vector<std::string> split (const std::string &s);

// Extract the second elements from a pair vector
pos_vector_t extract_second_from_pair(pos_pair_vector_t pos_pair_vector);

//hashing function for DnaString
unsigned dna2int(seqan::DnaString seq);



// Count the number of seeds in a mapping
unsigned count_seeds(read2pos_map_t results, unsigned current_read );

// Perform a Longest Increasing Sequence computation on a pair vector.
pos_pair_vector_t LIS_Pair(pos_pair_vector_t & pv);

// LIS for pair vector and spaced occurences
pos_pair_vector_t spaced_LIS_Pair(pos_pair_vector_t & pv, uint8_t k);

// Try to tell if a read may be a chimera. Return a score between 0 and 1.
float is_chimera(read2pos_map_t & mapping, unsigned read_len);


// Filtering of ressearch result before export
void filter_edges(read2pos_map_t & results, const seq_set_t & sequences, unsigned current_read, bool rev_comp, ef_params & params);