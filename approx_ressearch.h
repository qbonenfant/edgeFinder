#pragma once

// Include common types required for searching using SeqAn OSS
#include "types.h"
#include "params.h"
// Number of allowed errors
#define NB_ERR 1

// Main function for ressearch
void find_kmers(index_t & index, const seqan::DnaString & source_sequence, read2pos_map_t & pos_map, unsigned current_read, bool reversed, ef_params & params);