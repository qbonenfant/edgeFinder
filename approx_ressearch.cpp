#include <iostream>

//Include macros and required types
#include "approx_ressearch.h"
//#include "seedlib/seedlib.hpp"


using namespace seqan;

// Return the "complexity score" of a sequence
float lowComplexityScore(DnaString & sequence){
    // This approach is based on "DUST2" method,
    // but using 2-mers instead of 3-mers,
    // because we use far shorter sequences.
    unsigned l = length(sequence);
    unsigned counts[16] = { 0 };

    // reading using sliding window of 2
    for(unsigned i = 0; i < l-1; i++){
        uint8_t c = (uint8_t(sequence[i+1]) << 2) | uint8_t(sequence[i]);
        counts[c]++;
    }
    float score = 0;
    size_t sum = 0;
    for(auto v:counts){
        sum +=  v * (v-1);  
    }
    score =  sum / float(2 * (l-2));
    return score;
}


void find_kmers(index_t & index, const DnaString & source_sequence, read2pos_map_t & pos_map, unsigned current_read, bool reversed, ef_params & params){

	assert(params.ks >= 1);// avoid "infinite loop";
	
	// current position in the kmer list / sequence
	unsigned current_pos;

	// delegate function aggegating results
	 auto delegateParallel = [&](auto & iter, const DnaString & needle, int errors)
	        {
	        	// silencing warnings on unused template parameters
	        	(void)needle;
	        	(void)errors;
	            for (auto occ : getOccurrences(iter)){
	                // Identifying read Id and position on read
	                read_id_t  read_id  =  getValueI1(occ);
	                read_pos_t read_pos =  getValueI2(occ);
	                // checking we are not on the same read as source
		            if(current_read != read_id){    
		                // Checking if the match is not redundant.
	                    if( pos_map[read_id].empty() or pos_map[read_id].back().first != current_pos or 
	                    	std::abs(int64_t(pos_map[read_id].back().second) - read_pos) > NB_ERR + 1 ){
	                        pos_map[read_id].emplace_back(current_pos, read_pos);
	                    }
		            }        
				}	        
	        };

	// Going through sequence k-mers, jumping 'kmer_skip' position each time
	
	for(current_pos = 0; current_pos < length(source_sequence) - params.k + 1; current_pos += params.ks){
			DnaString kmer = infix(source_sequence, current_pos, current_pos + params.k);
			
			// Checking sequence complexity before searching
			if(lowComplexityScore(kmer) < params.lc){
				if(reversed){
					reverseComplement(kmer);
				}
				// Actual approximate search, filling the input position map
				find<0, NB_ERR>(delegateParallel, index, kmer , EditDistance() );
			}
	}
}