#include "utils.h"

#define DEBUG_FLAG 1

#if DEBUG_FLAG and not NDEBUG
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


//---------- DISPLAY AND EXPORT---------------

// Shortcut to print text in stdout with time stamp
const auto boot_time = std::chrono::steady_clock::now();
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}

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


//-------------- FASTA PARSING --------------
void parse_fasta(std::string filename, fasta_pair & fasta, uint8_t v){
    
    if(v>=1)
        print("PARSING FILE");

    SeqFileIn seqFileIn(toCString(filename));
    readRecords( fasta.first, fasta.second, seqFileIn);

    if(v>=2){
        print("Number of read parsed:");
        print(length(fasta.first));
    }

    if(v>=1)
        print("DONE");
}
//-------------- ############## --------------


//-------------- INDEX CREATION --------------
void create_index(std::string index_file, index_t & index, uint8_t v){

    bool success;
    if(index_file != ""){
        if(v>=1)
            print("LOADING FROM INDEX FILE");
        success = open(index,index_file.c_str());
        
    }
    else{
        if(v>=1)
            print("NO INDEX FILE, CREATING FROM SCRATCH...");

        success = indexCreate(index);
    }
    
    if(not success){
        print("INVALID INDEX FILE: EITHER PATH IS WRONG OR INPUT FILE CAN'T BE INDEXED");
        
    }

    else if(v>=1){
        print("DONE");
    }
    
}
//-------------- ############## --------------


// Extract the second elements from a pair vector
pos_vector_t extract_second_from_pair(pos_pair_vector_t pos_pair_vector){
    pos_vector_t second_part;
    second_part.reserve(pos_pair_vector.size());
    
    for(auto el : pos_pair_vector){
        second_part.emplace_back(el.second);
    }
    return(second_part);
}


// go through an associative array and build up the sum of the values, regardless of their key
template<typename TIterableNumbers>
inline unsigned vec_sum(TIterableNumbers number_list){
    unsigned total = 0;
    for(auto e: number_list){
        total += e;
    }
    return(total);
}


// LIS for vector of pairs. LIS is computed on second element
pos_pair_vector_t LIS_Pair(pos_pair_vector_t pv){
    pos_vector_t X = extract_second_from_pair(pv);
    const int N = X.size();
    int L = 0;
    pos_vector_t P(N);
    pos_vector_t M(N+1);
    

    for(int i=0; i<N; i++){
        int lo = 1;
        int hi = L;
        while(lo <= hi){
            int mid = ceil((lo+hi)/2);
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

// Try to tell if two reads are close enough to be isoforms or not
// This step could (should ?) be integrated directely during seeds processing,
// but is kept separated for modularity raisons.
// Return value is a small unsigned it with value between 0 and 3 (for now)
// 0 mean incoherent seed, 1 coherent seed but not enough cover, and 3 mean fine.
// 2 is used for the case where one read do have enough cover but not the second one
uint8_t is_iso(pos_pair_vector_t & pos_list,  unsigned l1, unsigned l2, uint8_t k, uint8_t ks, float max_diff_rate, float min_cover){


    std::vector<bool> coverage_ref(l1,0);
    std::vector<bool> coverage_tgt(l2,0);

    read_pos_t last_ref = pos_list[0].first;
    read_pos_t last_tgt = pos_list[0].second;

    // computing coverage and testing seed coherence
    int64_t dif_ref;
    int64_t dif_tgt;
    int64_t biggest_dif;
    int64_t smallest_dif;
    int64_t max_dif;
    int64_t relative_dif;

    for(auto p: pos_list){
        


        // coverage is not just on kmer start, but on all it's length
        for(int j=0; j< k; ++j){
            if(p.first + j < l1){
                coverage_ref[p.first + j] = true;
            }
            if(p.second + j < l2){
                coverage_tgt[p.second + j] = true;
            }
        }
        // how much did we advanced on reference and target ?
        dif_ref = std::abs(int64_t(p.first) - last_ref);
        dif_tgt = std::abs(int64_t(p.second) - last_tgt);

        // find the biggest gap
        biggest_dif = std::max(dif_ref, dif_tgt);
        // and smallest
        smallest_dif = std::min(dif_ref, dif_tgt);

        // maximum allowed difference
        max_dif= std::max(uint8_t(2), ks);
        if(smallest_dif != 0){
            max_dif = unsigned(std::max(int64_t(smallest_dif * max_diff_rate), max_dif));
        }

        // offset difference  between the two read
        relative_dif = std::abs(dif_ref - dif_tgt);

        // If we go over limit, just say it is not an isoform, for we can not be sure.
        if(biggest_dif > max_dif){
            DEBUG(biggest_dif);
            DEBUG(max_dif);
            DEBUG("INCOHERENT SEED");
            return(0);
        }

        // update the last positions
        last_ref = p.first;
        last_tgt = p.second;

    }
    // If seeds seem coherent, check if cover is high enough
    float reference_cover = vec_sum(coverage_ref) / float(l1);
    float target_cover    = vec_sum(coverage_tgt) / float(l2);
    
    DEBUG(std::to_string(reference_cover) + ", " + std::to_string(target_cover));

    if(reference_cover < min_cover and target_cover < min_cover){
        DEBUG("UNSUFFICIENT COVER ON BOTH READ");
        return(1);
    }
    else if(target_cover < min_cover and reference_cover >= min_cover){
        DEBUG("UNSUFFICIENT TARGET COVER");
        return(2);
    }       
    else if(target_cover >= min_cover and reference_cover < min_cover){
        DEBUG("UNSUFFICIENT REFERENCE COVER");
        return(2);
    }
        
    DEBUG("VALID ISOFORM");
    return(3);
}


// Converts DnaString to int. sequence longer than 32 will cause overflow.
inline unsigned dna2int(DnaString seq){
    
    unsigned value = 0;
    for(auto c : seq){
        value = value << 2 | uint8_t(c);    
    }
    return(value);
}

// adjust threshold value to kmer size
float adjust_threshold(float c_old, uint8_t k_old, uint8_t k_new ){
    float c_new = c_old * float( std::pow(k_new - 2 + 1,2) /  std::pow(k_old - 2 + 1,2));
    return(c_new);
}

// Just sum up the elements of an array
template<typename TIterable>
unsigned array_sum(TIterable int_array){
   
     unsigned sum = 0;
    for(auto el: int_array){
        sum+= el;
    }

    return(sum);
}



//-------------- SEED FILTERING ---------------

pos_pair_vector_t filter_seeds(pos_pair_vector_t & results, unsigned allocated_size, uint8_t ks){
    
    pos_pair_vector_t filtered_results;

    // Filtering redundant seeds.
    // allocating at least the same size as the original result array
    filtered_results.reserve(allocated_size);

    // keeping previous positions in memory
    unsigned last_first  = -1;
    unsigned last_second =  0;

    // for each seed occurences
    for(auto pos: results){
        // If it is a different seed, keep it 
        if( filtered_results.empty() or last_first != pos.first ){
            filtered_results.emplace_back(pos);
        }
        // else check that the second position do is different
        // by at least kmer skip size, minus the max error rate.
        else if(std::abs(int64_t(pos.second) - last_second ) > ks - 2 ){
            filtered_results.emplace_back(pos);
        }
        //updating last position
        last_first  = pos.first;
        last_second = pos.second;
    }

    return(filtered_results);
}