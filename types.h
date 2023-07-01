#pragma once

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <stdexcept>
#include <fstream>


// Application wide used types

// Fasta related types
using seq_id_set_t = seqan::StringSet<seqan::CharString>; //  Sequence id set type
using seq_set_t = seqan::StringSet<seqan::DnaString>; // Sequences set type
// Describe a fasta file using a pair of SeqAn data structures.
using fasta_pair_t = std::pair<seq_id_set_t, seq_set_t >; 

// Dataset info
using read_sizes_t = std::vector<unsigned>;


// Positions related types
using read_id_t = uint32_t;  // type for the id of the read. 16 bit were too short for files with >2M reads.
using read_pos_t = uint32_t; // type of the read position. Limited to 32 bits (position won't exceed 4 Gb)
using pos_vector_t = std::vector<read_pos_t>;// used to store occurences position
using read_pos_pair_t = std::pair<read_pos_t, read_pos_t>; // Position pair, used to describe mapping
using pos_pair_vector_t = std::vector<read_pos_pair_t>;    // Vector of position pairs
using bool_vector_t = std::vector<bool>;                   // Vector of boolean, for assignation tracking


// Index type, bidirectionnal FM-Index here.
using TFastConfig =  seqan::FastFMIndexConfig<void, size_t, 2, 0> ;
using index_t = seqan::Index<seq_set_t, seqan::BidirectionalIndex<seqan::FMIndex<void,TFastConfig> > >;

// Result storing type. Associate a read to a vector of read position pair.
using read2pos_map_t = std::unordered_map<read_id_t, pos_pair_vector_t >;





