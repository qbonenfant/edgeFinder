#pragma once

// required types
#include "types.h"
#include <seqan/seq_io.h>


//--------------------- DISPLAY----------------------------
// print template, containing a timer from program start
template<typename TPrintType>
void print(TPrintType text);


// Adjust the low complexity threshold to k-mer size, usin a k_old as base
float adjust_threshold(float c_old, uint8_t k_old, uint8_t k_new );

// Perform a Longest Increasing Sequence computation on a pair vector.
pos_pair_vector_t LIS_Pair(pos_pair_vector_t pv);

// Parse the input fasta file
void parse_fasta(std::string filename, fasta_pair & fasta, uint8_t v);

// Create the FM-Index
void create_index(std::string index_file, index_t & index, std::string index_folder , uint8_t v);

//hashing function for DnaString
unsigned dna2int(seqan::DnaString seq);
// LIS for pair vector
pos_pair_vector_t LIS_Pair(pos_pair_vector_t pv);

// check if two reads are from the same isoform
uint8_t is_iso(pos_pair_vector_t & pos_list,  unsigned l1, unsigned l2, uint8_t k, uint8_t ks, float max_diff_rate, float min_cover);

// Split strings on delimiter
std::vector<std::string> split (const std::string &s);