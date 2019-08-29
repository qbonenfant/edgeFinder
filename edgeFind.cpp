#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include <iostream>
#include <algorithm> 
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <unordered_set>


using namespace seqan;


// -------------------------- Defining common types ---------------------------------------------------------

// Index parameters, those are the default ones. // Maybe i should convert this typedef to using too...
typedef FastFMIndexConfig<void, size_t, 2, 0> TFastConfig;


// Fasta related types
using seq_id_set_t = StringSet<CharString>; //  Sequence id set type
using seq_set_t = StringSet<DnaString>; // Sequences set type
// Describe a fasta file using a pair of SeqAn data structures.
using fasta_pair = std::pair<seq_id_set_t, seq_set_t >; 

// Positions related types
using read_id_t = uint32_t;  // type for the id of the read. 16 bit were too short for files with >2M reads.
using read_pos_t = uint32_t; // type of the read position. Limited to 32 bits (position won't exceed 4 Gb)
using pos_vector_t = std::vector<read_pos_t>;// used to store occurences position
using read_pos_pair_t = std::pair<read_pos_t, read_pos_t>; // Position pair, used to describe mapping
using pos_pair_vector_t = std::vector<read_pos_pair_t>; //vector of position pairs


// Index type, bidirectionnal FM-Index here.
using index_t = Index<seq_set_t, BidirectionalIndex<FMIndex<void,TFastConfig> > >;

// Result storing type. Associate a read to a vector of read position pair.
using read2pos_map_t = std::map<read_id_t, pos_pair_vector_t >;


// ----------------------------------------------------------------------------------------------------------


//-------------- UTILS FUNCTION --------------

// Shortcut to print text in stdout with time stamp
const auto boot_time = std::chrono::steady_clock::now();
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
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

// go through an associative array and build up the sum of the values, regardless of their key
template<typename TIterableNumbers>
inline unsigned map_sum(TIterableNumbers number_list){
    unsigned total;
    for(auto e: number_list){
        total += e.second;
    }
    return(total);
}

// Extract the first elements from a pair vector
pos_vector_t extract_first_from_pair(pos_pair_vector_t pos_pair_vector){
    pos_vector_t second_part;
    second_part.reserve(pos_pair_vector.size());
    
    for(auto el : pos_pair_vector){
        second_part.emplace_back(el.first);
    }
    return(second_part);
}

// return true if the sequence is considered "low complexity"
inline  bool haveLowComplexity(DnaString & sequence, float threshold){
    // New version, using "DUST2" method
    // scanning 2-mers, squaring count, discard if over limit
    unsigned l = length(sequence);

    unsigned counts[16] = { 0 };
    // reading using sliding window of 2
    for(int i = 0; i < l-1; i++){
        uint8_t c = (uint8_t(sequence[i+1]) << 2) | uint8_t(sequence[i]);
        counts[c]++;
        
    }
    float s = 0;
    size_t sum = 0;
    for(auto v:counts){
        sum +=  v * (v-1);  
    }
    s =  sum / float(2 * (l-2));
    return s>= threshold;
        
}


// Wikipedia implementation of an efficient LIS
// Works in O(n log n) time.
pos_vector_t LIS(pos_vector_t X){
    
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
    pos_vector_t S(L);
    int k = M[L];
    for(int i= L-1; i>=0; i--){
        S[i] = X[k];
        k = P[k];
    }
    return(S);
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
bool is_iso(pos_pair_vector_t & pos_list,  unsigned l1, unsigned l2, uint8_t k=16, uint8_t ks=3, float max_diff_rate=1.30, float min_cover= 0.75){

    std::map<read_pos_t, uint8_t> coverage_ref;
    std::map<read_pos_t, uint8_t> coverage_tgt;

    read_pos_t last_ref = pos_list[0].first;
    read_pos_t last_tgt = pos_list[0].second;

    for(auto p: pos_list){
        
        // computing coverage and testing seed coherence
        unsigned dif_ref;
        unsigned dif_tgt;
        unsigned biggest_dif;
        unsigned smallest_dif;
        unsigned max_dif;
        unsigned relative_dif;

        // coverage is not just on kmer start, but on all it's length
        for(int j=0; j< k; ++j){
            coverage_ref[p.first + j] = 1;
            coverage_tgt[p.second + j] = 1;
        }
        // how much did we advanced on reference and target ?
        dif_ref = std::abs(p.first - last_ref);
        dif_tgt = std::abs(p.second - last_tgt);

        // find the biggest gap
        biggest_dif = std::max(dif_ref, dif_tgt);
        // and smallest
        smallest_dif = std::min(dif_ref, dif_tgt);

        // maximum allowed difference
        max_dif= ks;
        if(smallest_dif != 0){
            max_dif = unsigned(std::max(unsigned(smallest_dif * max_diff_rate), max_dif));
        }

        // offset difference  between the two read
        relative_dif = std::abs(dif_ref - dif_tgt);

        // If we go over limit, just say it is not an isoform, for we can not be sure.
        if(biggest_dif > max_dif){
            return(false);
        }
    }
    // If seeds seem coherent, check if cover is high enough
    float reference_cover = map_sum(coverage_ref) / l1;
    float target_cover    = map_sum(coverage_tgt) / l2;
    if(reference_cover < min_cover or target_cover < min_cover){
        return(false);
    }
  
    return(true);
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
        else if(abs(pos.second - last_second ) > ks - 2 ){
            filtered_results.emplace_back(pos);
        }
        //updating last position
        last_first  = pos.first;
        last_second = pos.second;
    }

    return(filtered_results);
}


//-------------- EXPORT FUNCTION --------------

void export_read_result( read2pos_map_t results, read_id_t current_read_id, unsigned nk, bool reverse, const seq_id_set_t & realIds,  const seq_set_t & sequences,  std::ofstream & output_file){
    
    for( auto it=results.begin(); it != results.end(); ++it){
        
        
        // Avoid printing the same id, and results that are too short
        if(it->first != current_read_id and it->second.size() >= nk ){
            output_file << realIds[current_read_id];
            output_file << "\t" << length(sequences[current_read_id]);
            output_file << "\t" << realIds[it->first];
            output_file << "\t" << length(sequences[it->first]);
            output_file << "\t" << reverse;

            bool isoforms = is_iso(it->second, length(sequences[current_read_id]), length(sequences[it->first]), 1.25, 0.75);
            output_file << "\t" << isoforms;

            // Exporting seeds.
            for(auto pos: it->second){
                    output_file << "\t" << pos.first << "," << pos.second;    
            }
            output_file << "\n";
        }
    }
}



//-------------- ############## --------------


//-------------- FASTA PARSING --------------
void parse_fasta(std::string filename, fasta_pair & fasta, uint8_t v){
    
    if(v>=1)
        print("PARSING FILE");

    SeqFileIn seqFileIn(toCString(filename));
    readRecords( fasta.first, fasta.second, seqFileIn);

    if(v>=1)
        print("DONE");
}
//-------------- ############## --------------


//-------------- INDEX CREATION --------------
void create_index(std::string index_file, index_t & index, uint8_t v){

    if(index_file != ""){
        if(v>=1)
            print("LOADING FROM INDEX FILE");
        open(index,index_file.c_str());
        
    }
    else{
        if(v>=1)
            print("NO INDEX FILE, CREATING FROM SCRATCH...");
        indexCreate(index);
    }
    
    if(v>=1)
        print("DONE");
}
//-------------- ############## --------------


//--------------APPROX RESSEARCH--------------
void approxCount(const seq_id_set_t & ids, const seq_set_t  & sequences, index_t & index, std::ofstream & output_file, uint8_t k, uint8_t ks, uint8_t nk ,const uint8_t nb_thread, double lc, bool rc, bool sampling, uint8_t v){

    const unsigned NB_ERR = 1;

    if(v>=1)
        print("STARTING RESSEARCH");
   
    // Set of already processed reads. Needs to be shared between all thread.
    std::unordered_set<read_id_t> processed_reads;
    
    // Setting number of thread
    omp_set_num_threads(nb_thread);

    // Some counters for run statistics, one per thread to avoid colission.
    // A 64 bit unsigned int should be large enough to keep track of those.
    std::vector<unsigned>  nb_ressearch(nb_thread);
    std::vector<unsigned>  collisions(nb_thread);
    std::vector<unsigned>  total_kmer(nb_thread);
    std::vector<unsigned>  lc_kmer(nb_thread);
            

    // ######################################################################################################
    // PARALLEL SECTION #####################################################################################
    // ######################################################################################################
    #pragma omp parallel  
    {

        // INITIALISING VARIABLES FOR EACH THREAD

        // Mapping processing variables
        int direction = 0 ;         // direct or reverse direction
        read2pos_map_t direct_pos;  // store  pair of position for direct 
        read2pos_map_t reverse_pos; // and reverse complement kmers.
        
        
        
        // Keeping track of processed kmer (using hash)    
        std::unordered_map<unsigned, bool> processedKmer;
        
        // current position in the read, need to be declared here to be used inside the lambda.
        unsigned current_pos;
        // number of read
        unsigned nb_read = length(sequences);


        // Lambda used to process occurences. 
        // This function no longer use pragma omp critical since all variables here are 
        // supposed to be private.
        // Only the result export should be protected.
        auto delegateParallel = [&](auto & iter, const DnaString & needle, int errors)
        {
            for (auto occ : getOccurrences(iter)){
                // Identifying read Id and position on read
                read_id_t  readId  =  getValueI1(occ);
                read_pos_t readPos =  getValueI2(occ);
                

                // Checking if the match is not redundant.
                if(direction == 0){
                    if( direct_pos[readId].empty() or direct_pos[readId].back().first != current_pos ){
                        direct_pos[readId].emplace_back(current_pos, readPos);
                    }
                    else if( abs(direct_pos[readId].back().second - readPos) > NB_ERR + 1 ){
                        direct_pos[readId].emplace_back(current_pos, readPos);
                    }
                }

                else{
                    if( reverse_pos[readId].empty() or reverse_pos[readId].back().first != current_pos ){
                        reverse_pos[readId].emplace_back(current_pos, readPos);
                    }
                    else if( abs(reverse_pos[readId].back().second - readPos) > NB_ERR + 1  ){
                        reverse_pos[readId].emplace_back(current_pos, readPos);
                    }
                    
                }
            }        
        
        };

        // ##################################################################################################
        // MAIN LOOP ########################################################################################
        // ##################################################################################################
    
        // Going through the read list using a dynamic range allocation strategy
        #pragma omp for schedule(dynamic)
        for(read_id_t r=0; r<nb_read; r++)
        {
            // Progress tracking, if file is large enough.
            if(v>=2 and nb_read>=100 and (r) %(nb_read/100) == 0){
                print( std::to_string(round(float(r+1)/nb_read*100)) + "%" );
                print( r );
            }

            // Checking read id is not already found, if sampling is activated
            if( not sampling  or processed_reads.find(r) == processed_reads.end() ){
                
                // #pragma omp critical
                // {
                //     processed_reads;
                // }

                // cleaning result set
                direct_pos.clear();
                if(rc){
                    reverse_pos.clear();
                }
                //processedKmer.clear();
                DnaString readSequence = sequences[r];

                // ##########################################################################################
                // RESSEARCH ################################################################################
                // ##########################################################################################

                // going through the read using k size window, every ks
                for( current_pos = 0; current_pos < length(readSequence) - k + 1; current_pos += ks ){
                    
                    DnaString km = infix(readSequence, current_pos, current_pos + k);
                    // Counting redundant kmer search. 
                    // hashing kmer to find out if we searched it before.
                    unsigned kmHash = dna2int(km);
                    if(not processedKmer[kmHash]){
                        processedKmer[kmHash] = true;    
                    }
                    else{
                        collisions[omp_get_thread_num()] += 1;
                    }

                    // Counting the number of kmer we looked at
                    total_kmer[omp_get_thread_num()] += 1;

                    // Checking if the kmer has low complexity 
                    bool haveLc = haveLowComplexity(km,lc);

                    // If not, search it in the index.
                    if(not haveLc){

                        // Forward strand ressearch
                        direction = 0;
                        find<0, NB_ERR>(delegateParallel, index, km , EditDistance() );
                        
                        // Counting number of ressearch done
                        nb_ressearch[omp_get_thread_num()] +=1;

                        // If we want to also search for reverse complement
                        if(rc){
                            // ressearching using reverse complement of the k-mer too
                            reverseComplement(km);
                            direction = 1;
                            find<0, NB_ERR>(delegateParallel, index, km , EditDistance() );
                            
                            // Also counting ressearch in reverse
                            nb_ressearch[omp_get_thread_num()] +=1;
                        }
                    }
                    else{
                        // Counting low complexity kmers
                         lc_kmer[omp_get_thread_num()] +=1;
                    }
                    
                }

                // ##########################################################################################
                // FILTERING ################################################################################
                // ##########################################################################################

                // Checking number of matches in forward

                for(auto it = direct_pos.begin(); it != direct_pos.end(); ++it ){
                    // looking for the longet streak of kmers that are colinear to our read.
                    it->second = LIS_Pair(it->second);
                    it->second = filter_seeds(it->second, it->second.size(), ks);
                    if( it->second.size() >= nk  and is_iso(it->second, length(sequences[r]), length(sequences[it->first]), 1.25, 0.75)){
                        // if read id is associated, do not use it again
                        #pragma omp critical
                        processed_reads.insert(it->first);
                    }
                }

                if(rc){
                    for(auto it = reverse_pos.begin(); it != reverse_pos.end(); ++it ){
                        // Reversing the results, otherwise LIS won't really make sens.
                        std::reverse( it->second.begin(), it->second.end());
                        // Looking for the longet streak of kmers that are colinear to our read.
                        it->second = LIS_Pair(it->second);
                        // Filtering redundant seeds
                        it->second = filter_seeds(it->second, it->second.size(), ks);

                        if( it->second.size() >= nk  and is_iso(it->second, length(sequences[r]), length(sequences[it->first]), 1.25, 0.75)) {
                            // if read id is associated, do not use it again
                            #pragma omp critical
                            processed_reads.insert(it->first);
                        }
                    }
                }

                // EXPORT THE RESULTS
                #pragma omp critical
                {
                    // export results :  map, read id, threshold, reverse?, id list, output stream
                    export_read_result(direct_pos,  r, nk, false, ids, sequences, output_file );
                    export_read_result(reverse_pos, r, nk, true,  ids, sequences, output_file );
                }
                
            }
            // ##############################################################################################
            // END OF RESSEARCH #############################################################################
            // ##############################################################################################
        }
        // ##################################################################################################
        // END OF MAIN LOOP #################################################################################
        // ##################################################################################################
    }
    // ######################################################################################################
    // END OF PARALLEL SECTION ##############################################################################
    // ######################################################################################################


    if(v>=1){
        print("DONE");
    }

    if(v>=2){
        std::cout << "Total number of kmer searched: " <<  array_sum(total_kmer) << std::endl;
        std::cout << "Number of collisions: " << array_sum(collisions);
        std::cout << " ("  << array_sum(collisions)*100/array_sum(total_kmer) << "%)"     << std::endl;
        std::cout << "Number of discarded kmer (lc): " << array_sum(lc_kmer);
        std::cout << " ("  << array_sum(lc_kmer)*100/array_sum(total_kmer) << "%)"     << std::endl;
        std::cout << "Number of researchs: " << array_sum(nb_ressearch)  << std::endl;
        
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
        "i", "index", "path and prefix of the index files",
        seqan::ArgParseArgument::STRING, "index file"));

    addOption(parser, seqan::ArgParseOption(
        "o", "outFile", "path to the output file",
        seqan::ArgParseArgument::STRING, "output file"));


    // Parameters for the algorithm
    addSection(parser, "Ressearch parameters");    

    addOption(parser, seqan::ArgParseOption(
        "k", "k", "Size of the kmers",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "nk", "nbKmer", "Minimum number of common k-mer needed to create an edge between two reads.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "ks", "kmer-skip", "Limit ressearch to 1/ks k-mers in the read",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "lc", "low-complexity", "Threshold for low complexity",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, seqan::ArgParseOption(
        "rc", "revComp", "Ressearsh k-mer's rev-comp too. Set to 0 to deactivate.",
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
    std::string index_file = "";  // index prefix, if any                              
    unsigned nb_thread = 4;       // default number of thread (4)                      
    //unsigned nbErr = 2;         // default number of errors EDIT: CAN NOT BE ASSIGNED
    unsigned v = 1;               // verobisty, default = 1                            
    unsigned k = 16;              // kmerSize, default = 16                            
    unsigned ks = 3;              // k-mer skip size, default = 3                      
    unsigned nk= 10;              // Minimum common k-mer, default = 10                 
    double   lc= 1.25;            // "dust2" low complexity threshold, default = 1.25  
    bool     rc =  true;          // Rev Comp research ? True by default               
    bool sampling =  true;        // Sampling method ? True by default                 
    // --------------------------------------------------------------------------------



    // ASSIGNING VALUES FROM PARSER ---------------------------------------------------
    //
    // Boolean values
    // checking if revcomp is activated
    if( isSet(parser, "revComp") ){
        int val;
        getOptionValue(val, parser, "rc");
        rc = val!=0;
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
    getOptionValue(index_file, parser, "i");
    
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
    output_file << " index:" << index_file;
    output_file <<std::endl;
    output_file <<  "output:" << output << std::endl;

    // PROGRAM STARTING POINT

    // Fasta parsing
    fasta_pair fasta;
    parse_fasta(fasta_file, fasta, v);

    // Index creation
    
    if(v>=1)
        print("INITIALISING INDEX");
    index_t index(fasta.second);
    create_index(index_file, index, v);

    // Ressearch and export
    approxCount(fasta.first, fasta.second, index, output_file, k, ks, nk, nb_thread, lc, rc, sampling, v);
    
    
    // closing

    output_file.close();
    if(v>=1){
        print("PROGRAM END");
    }

    return 0;
}
//-------------- ##### END ##### --------------