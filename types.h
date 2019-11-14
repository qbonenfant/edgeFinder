#pragma once

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>


// Application wide used types

// Fasta related types
using seq_id_set_t = seqan::StringSet<seqan::CharString>; //  Sequence id set type
using seq_set_t = seqan::StringSet<seqan::DnaString>; // Sequences set type
// Describe a fasta file using a pair of SeqAn data structures.
using fasta_pair = std::pair<seq_id_set_t, seq_set_t >; 

// Positions related types
using read_id_t = uint32_t;  // type for the id of the read. 16 bit were too short for files with >2M reads.
using read_pos_t = uint32_t; // type of the read position. Limited to 32 bits (position won't exceed 4 Gb)
using pos_vector_t = std::vector<read_pos_t>;// used to store occurences position
using read_pos_pair_t = std::pair<read_pos_t, read_pos_t>; // Position pair, used to describe mapping
using pos_pair_vector_t = std::vector<read_pos_pair_t>; //vector of position pairs


// Index type, bidirectionnal FM-Index here.
using TFastConfig =  seqan::FastFMIndexConfig<void, size_t, 2, 0> ;
using index_t = seqan::Index<seq_set_t, seqan::BidirectionalIndex<seqan::FMIndex<void,TFastConfig> > >;

// Result storing type. Associate a read to a vector of read position pair.
using read2pos_map_t = std::map<read_id_t, pos_pair_vector_t >;

// Keeping track of isoforms status between target and source
using isoform_map_t = std::map<unsigned, uint8_t>;

// storing node type: repr or delegate.
using node_type_t = std::map<std::string, bool>;

// storing node degree
using node_degree_t = std::map<std::string, unsigned>;