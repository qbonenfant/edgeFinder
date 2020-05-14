#include "utils.h"

//#define EF_DEBUG 1

#if EF_DEBUG and not NDEBUG
#define DEBUG(x) std::cout << x << "\n"
#else
#define DEBUG(x)
#endif

using namespace seqan;

//-------------- UTILS FUNCTION --------------
//  This file contains miscellaneous function
//  used for seeds processing, conversions
//  index creation, file loading and export.
//--------------------------------------------



///////////////////////////////////////////////////////////////////////////////
// ------------------------  DISPLAY AND EXPORT -----------------------------//

// Shortcut to print text in stdout with time stamp
const auto boot_time = std::chrono::steady_clock::now();
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}

// Export the mapping data to the outfile in .edges format
void export_edges(read2pos_map_t & results, read_id_t current_read, const seq_id_set_t & realIds, const seq_set_t & sequences, bool reverse, float score, std::ofstream & output_file){
    
    for(auto it = results.begin(); it != results.end(); ++it){

        output_file << realIds[current_read];
        output_file << "\t" << length(sequences[current_read]);
        output_file << "\t" << realIds[it->first];
        output_file << "\t" << length(sequences[it->first]);
        output_file << "\t" << reverse;
        output_file << "\t" << score;

        // Exporting seeds.
        for(auto pos: it->second){
                output_file << "\t" << pos.first << "," << pos.second;
        }
        output_file << "\n";
    }
}

//---------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
//---------------------------- FASTA PARSING --------------------------------//

// Parse a fasta file using SeqAn structure.
// Data container is passed by reference and filled.
void parse_fasta(fasta_pair_t & fasta, ef_params & params){
    
    seqan::SeqFileIn seqFileIn(seqan::toCString(params.fasta_file));
    readRecords( fasta.first, fasta.second, seqFileIn);
}
//---------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
//----------------------------- INDEX CREATION ------------------------------//

void create_index(index_t & index, ef_params & params){
    // TODO: Fix the output message and fail management.
    // Errors should go to stderr, and index creation
    // can not fail, or else  the program needs to end.

    bool success = false;
    bool save_success = false;
    if(params.index_file != ""){
        if(params.v>=1)
            print("LOADING FROM INDEX FILE");
        success = open(index, params.index_file.c_str());
        if(params.index_folder != ""){
            print("LOADED INDEX WON'T BE SAVED, IT ALREADY EXISTS");
        }
        
    }
    else{
        if(params.v>=1)
            print("NO INDEX FILE, CREATING FROM SCRATCH...");

        success = indexCreate(index);
    }
    
    if(not success){
        print("INVALID INDEX FILE: EITHER PATH IS WRONG OR INPUT FILE CAN'T BE INDEXED");
        
    }

    if(params.index_file == "" and params.index_folder != ""){
        if(params.v>=1){
            print("SAVING CREATED INDEX IN SPECIFIED LOCATION.");
        }
        save_success = save(index, seqan::toCString(params.index_folder));
        
        if(save_success){
            print("INDEX SAVED!");
        }
        else if (params.v>=1){
            print("SOMETHING HAS GONE WRONG WHILE SAVING INDEX, CONTINUING WITHOUT SAVING.");
        }
    }
    print("DONE");
}
// --------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// --------------------------- MISC FUNCTIONS -------------------------------//

// Spliting strings on delimiter
std::vector<std::string> split (const std::string & s) {
    char delim = *"\t";
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }
    return result;
}

// Extract the second elements from a pair vector
pos_vector_t extract_second_from_pair(pos_pair_vector_t pos_pair_vector){
    pos_vector_t second_part;
    second_part.reserve(pos_pair_vector.size());
    
    for(auto el : pos_pair_vector){
        second_part.emplace_back(el.second);
    }
    return(second_part);
}

// Go through an iterable and build up the sum of the values
template<typename TIterableNumbers>
unsigned iter_sum(TIterableNumbers number_list){
    unsigned total = 0;
    for(auto e: number_list){
        total += e;
    }
    return(total);
}

// Converts DnaString to int. sequence longer than 32 will cause overflow.
unsigned dna2int(seqan::DnaString seq){

    unsigned value = 0;
    if(length(seq)>32){
        throw std::invalid_argument("K-mer can not exceed 32 bases for uint64_t conversion.");
    }

    for(auto c : seq){
        value = value << 2 | uint8_t(c);    
    }
    return(value);
}



// Count the number of seeds in a mapping
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

// --------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
//------------------- LONGEST INCREASING SUBSEQUENCE-------------------------//

// LIS for vector of pairs. LIS is computed on second element
pos_pair_vector_t LIS_Pair(pos_pair_vector_t & pv){
    pos_vector_t X = extract_second_from_pair(pv);
    const int N = X.size();
    int L = 0;
    pos_vector_t P(N);
    pos_vector_t M(N+1);
    

    for(int i=0; i<N; i++){
        int lo = 1;
        int hi = L;
        while(lo <= hi){
            int mid = (lo+hi+1)/2;
            if(X[M[mid]] <= X[i])
                lo = mid+1;
            else
                hi = mid-1;
        }
        int newL = lo;
        P[i] = M[newL-1];
        M[newL] = i;

        if(newL > L)
            L = newL;
    }
    pos_pair_vector_t S(L);
    int k = M[L];
    for(int i= L-1; i>=0; i--){
        S[i] = pv[k];
        k = P[k];
    }
    return(S);
}


// Same as LIS pair, but element need to be spaced by km_size
pos_pair_vector_t spaced_LIS_Pair(pos_pair_vector_t & pv, uint8_t km_size){    

    pos_vector_t X = extract_second_from_pair(pv);
    const int N = X.size();
    int L = 0;
    pos_vector_t P(N);
    std::vector<int> M(N+1, -1);
    M[0] = 0;    

    for(int i=0; i<N; i++){
        int lo = 1;
        int hi = L;
        while(lo <= hi){
            int mid = (lo+hi+1)/2;
            if(X[M[mid]] + km_size <= X[i])
                lo = mid+1;
            else
                hi = mid-1;
        }
        int newL = lo;
        P[i] = M[newL-1];
        
        // is there a value for M[newL] ?
        if(M[newL] >= 0){
            // is value for newL really better ?
            if(X[M[newL]] >= X[i] ){
                // is it compatible with interval?
                if( newL==1 or X[P[i]] + km_size <= X[i] ){
                    M[newL] = i;
                }
            }
        }
        // If no value, just set it;
        else{
            M[newL] = i;
        }

        if(newL > L)
            L = newL;
    }
    pos_pair_vector_t S(L);
    int k = M[L];
    for(int i= L-1; i>=0; i--){
        S[i] = pv[k];
        k = P[k];
    }
    return(S);
}
// --------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// ------------------------ CHIMERA DETECTION -------------------------------//

// Try to tell if a read may be a chimera. Return a score between 0 and 1.
float is_chimera(read2pos_map_t & mapping, unsigned read_len){

    // mapping seed distribution for each read
    std::map<read_pos_t, std::map<read_id_t, float> > AB_map;
    // curent AB value for each read. Starts at 1.
    std::unordered_map<read_id_t, float> current_read_val;
    // number of read in the mapping
    unsigned nb_read = mapping.size();

    // Creating lookup structures for faster computation of the score
    for(auto read_seeds: mapping){
        // init the starting AB value for each read
        current_read_val[read_seeds.first] = 1.0;
        // number of seeds on this read
        int r_nb_seed = read_seeds.second.size();
        int current_seed_nb = 1;
        // for each seed position, store the AB diff
        for(read_pos_pair_t pos_pair: read_seeds.second){
            float AB_score = abs((2 * current_seed_nb) - r_nb_seed ) / r_nb_seed;
            AB_map[pos_pair.first][read_seeds.first] = AB_score; 
            current_seed_nb ++;
        }
    }

    // Computing score for each seed position and keeping max
    float max_score = 0;
    for(auto s_pos: AB_map){
        
        // expected AB difference at this seed position for the source read
        float exp = std::abs((2 * s_pos.first) - read_len) / read_len;

        // current score at a position is the sum of the diff between AB score and exp
        // normalized by the number of read mapped
        float current_score = 0.0;

        // for each mapped read, fetch the AB score
        for(auto rid: current_read_val){
            // checking if read id has a value in the s_pos map
            if(s_pos.second.find( rid.first ) != s_pos.second.end()){
                // If yes, set current AB value for this read
                rid.second = s_pos.second[rid.first];
            }
            current_score += rid.second - exp;
        }
        current_score /= nb_read;
        max_score = std::max(max_score, current_score);
    }
    return(max_score);
}
// --------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// -------------------------- SEEDS PRE-FILTER ------------------------------//

// DEPRECATED: This functionnality is now included in
// the delegate function of the ressearch section

// This function ensure that we do not have duplicate seeds before starting the LIS
// pos_pair_vector_t filter_seeds(pos_pair_vector_t & pair_vector, unsigned allocated_size, uint8_t ks){
    
//     pos_pair_vector_t filtered_pair_vector;

//     // Filtering redundant seeds.
//     // allocating at least the same size as the original result array
//     filtered_pair_vector.reserve(allocated_size);

//     // keeping previous positions in memory
//     unsigned last_first  = -1;
//     unsigned last_second =  0;

//     // for each seed occurences
//     for(auto pos: pair_vector){
//         // If it is a different seed, keep it 
//         if( filtered_pair_vector.empty() or last_first != pos.first ){
//             filtered_pair_vector.emplace_back(pos);
//         }
//         // else check that the second position do is different
//         // by at least kmer skip size, minus the max error rate.
//         else if(std::abs(int64_t(pos.second) - last_second ) > ks - 2 ){
//             filtered_pair_vector.emplace_back(pos);
//         }
//         //updating last position
//         last_first  = pos.first;
//         last_second = pos.second;
//     }
//     return(filtered_pair_vector);
// }
// --------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// ---------------------------- EDGES  FILTER -------------------------------//
void filter_edges(read2pos_map_t & results, const seq_set_t & sequences, unsigned current_read, bool rev_comp, ef_params & params){
    
    // manually moving the iterator makes for easier map modification.
    auto result_iterator = results.begin();

    // Checking each mapped read
    while(result_iterator != results.end()){

        // If we work with opposite orientation reads reverse seed vector first
        if(rev_comp){
            std::reverse( result_iterator->second.begin(), result_iterator->second.end());
        }

        // Keeping longest colinear streak of seeds with no overlap
        if(params.lis_mode == 1 ){
            result_iterator->second = spaced_LIS_Pair(result_iterator->second, params.k);
        }
        // Or with overlap allowed, depending on option
        else{
            result_iterator->second = LIS_Pair(result_iterator->second);
        }

        // fetching reads length and deriving required number of common k-mers
        unsigned l1 = length(sequences[current_read]);
        unsigned l2 = length(sequences[result_iterator->first]);
        unsigned rnk = std::max(params.nk, (unsigned)std::floor(params.kp * std::min(l1,l2) / params.k )  );

        // If we do not have enough seeds, reject the mapping
        if(result_iterator->second.size() < rnk){
            result_iterator = results.erase(result_iterator);
        }
        else{
            // else just fetch the next element
            result_iterator ++;
        }       
    }
}
// --------------------------------------------------------------------------//
///////////////////////////////////////////////////////////////////////////////
