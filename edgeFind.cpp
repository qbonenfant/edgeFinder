#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include <iostream>
#include <cstdlib>
#include <fstream>


using namespace seqan;


// Setting the index
typedef FastFMIndexConfig<void, size_t, 2, 0> TFastConfig;

using read_pos_t = uint16_t;
using read_id_t = unsigned;
using pos_vector_t = std::vector<read_pos_t>;
using read2pos_map_t = std::map<int,std::vector<read_pos_t> >;
using index_t = Index<StringSet<DnaString>, BidirectionalIndex<FMIndex<void,TFastConfig> > >;
using int_vector = std::vector<uint16_t>;

const auto boot_time = std::chrono::steady_clock::now();
// Shortcut to print text in stdout
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}


inline unsigned dna2int(DnaString seq){
    
    unsigned value = 0;
    for(auto c : seq){
        value = value << 2 | uint8_t(c);    
    }
    return(value);
}


float adjust_threshold(float c_old, uint8_t k_old, uint8_t k_new ){
    float c_new = c_old * float((k_new - 2 + 1)^2 / (k_old - 2 + 1)^2);
    return(c_new);
}

int_vector LIS(int_vector X){
    // Wikipedia implementation of an efficient LIS
    // Works in O(n log n) time.
    // Bench test shows it can compute LIS for 10^9 sequences of 5000 elements in less than a second.
    // Good enough for me.
    const int N = X.size();
    int L = 0;
    int_vector P(N);
    int_vector M(N+1);  
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
    int_vector S(L);
    int k = M[L];
    for(int i= L-1; i>=0; i--){
        S[i] = X[k];
        k = P[k];
    }
    return(S);
}

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



// Function used to print the "clustering" results
void printReadNames(int currentReadId, read2pos_map_t& readIdMap , StringSet<CharString>& realIds ,std::ofstream &outputFile)
{
    
    outputFile << realIds[currentReadId];
    for(auto it=readIdMap.begin(); it != readIdMap.end(); ++it){

        int id =  it->first;
        // Avoid printing the same id
        if(id != currentReadId){
            outputFile << "\t" << realIds[id];
            outputFile << "\t" << it->second.size();
        }
    }

    outputFile << "\n";
}

// Function used to print the read names associated to current read, if the length of the LIS of mapped pos is high enough 
void printReadNamesWithLIS(int currentReadId, read2pos_map_t& readIdMap , StringSet<CharString>& realIds ,std::ofstream &outputFile)
{
    
    outputFile << realIds[currentReadId];
    for(auto it=readIdMap.begin(); it != readIdMap.end(); ++it){

        int id =  it->first;
        // Avoid printing the same id
        if(id != currentReadId){
            outputFile << "\t" << realIds[id];
            outputFile << "\t" << LIS(it->second).size();
        }
    }

    outputFile << "\n";
}


// Function used to print all the positions for associated reads
void printReadNamesAndPos(int currentReadId, read2pos_map_t& readIdMap , StringSet<CharString>& realIds ,std::ofstream &outputFile)
{
    
    for(auto it=readIdMap.begin(); it != readIdMap.end(); ++it){
        int id =  it->first;
        // Avoid printing the same id
        if(id != currentReadId){
            outputFile << realIds[currentReadId];
            outputFile << "\t" << realIds[id];
            for(auto pos: it->second){
                outputFile << "\t" << pos;
            }
            outputFile << "\n";
        }
    }

    outputFile << "\n";
}


// Search and count a kmer list in a fasta file, at at most a levenstein distance of 2.
void approxCount(const std::string& filename, const std::string & indexFile, const int k, const int& nb_thread, std::ofstream &outputFile, int v, int bloc_size, int treshold, bool rc, double lc, bool sampling ){
    
    // Max number of errors, need to be fixed at compile time
    const uint8_t NB_ERR = 1;

    if(v>=1)
        print("PARSING FILE");

    // Parsing input fasta file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn(toCString(filename));
    readRecords(ids, seqs, seqFileIn);

    if(v>=1)
        print("DONE");

    // Constants (should be in upper case)
    const int nbRead  = length(ids);     // Number of reads in the file

    
    // INDEX CREATION ----------
    if(v>=1)
        print("CREATING INDEX");

    index_t index(seqs);

    if(indexFile != ""){
        if(v>=1)
            print("LOADING FROM INDEX FILE");
        open(index,indexFile.c_str());
        
    }
    else{
        if(v>=1)
            print("NO INDEX FILE, CREATING FROM SCRATCH...");
        indexCreate(index);
    }
    
    if(v>=1)
        print("DONE");
    
    // -------------------------
    // RESSEARCH ---------------
    // -------------------------

    if(v>=1)
        print("STARTING RESSEARCH");
   
    std::set<read_id_t> stopSearch;
    //unsigned searchCount[nb_thread]={0};
    
    omp_set_num_threads(nb_thread);
    #pragma omp parallel  shared(ids, seqs, treshold, rc, stopSearch, lc, sampling, index)
    {
        // storing results
        std::array<read2pos_map_t, 2> results;
        int direction = 0 ;
        read2pos_map_t keptId;
        std::unordered_map<unsigned, bool> processedKmer;

        

        auto delegateParallel = [&](auto & iter, const DnaString & needle, int errors)// __attribute__((noinline))
        {
            for (auto occ : getOccurrences(iter)){
                // Identifying read Id and position on read    
                // and testing if the occurence is on one read only
                // (not between two reads)
                int readId= getValueI1(occ);
                read_pos_t readPos =  getValueI2(occ); 
                #pragma omp critical
                if(std::find(results[direction][readId].begin(), results[direction][readId].end(), readPos ) == results[direction][readId].end() ){
                   results[direction][readId].push_back(readPos);
                }
            }
        };

        // going through the read list
        #pragma omp for schedule(dynamic)
        for(read_id_t r=0; r<nbRead; r++)
        {

            if(v>=2 and nbRead>=100 and  (r+1) %(nbRead/100) == 0){

                print( std::to_string(round(float(r+1)/nbRead*100)) + "%" );
                print( r );
            }

            // checking read id is not already found, if sampling is activated
                    
            if( not sampling  or stopSearch.find(r) == stopSearch.end()){

                // cleaning result set
                results[0].clear();
                results[1].clear();
                processedKmer.clear();
                DnaString readSequence = seqs[r];
                // going through the read using K size window
                for( int i = 0; i < length(readSequence) - k + 1; i+= bloc_size ){
                    
                    DnaString km = infix(readSequence, i, i+k);
                    
                    // hashing kmer to find out if we searched it before.
                    unsigned kmHash = dna2int(km);
                    if(not processedKmer[kmHash]){
                        processedKmer[kmHash] = true;
                        bool haveLc = haveLowComplexity(km,lc);
                        if(not haveLc){
                            // Forward strand ressearch
                            direction = 0;
                            find<0, NB_ERR>(delegateParallel, index, km , EditDistance() );
                            
                            // Reverse complement
                            if(rc){
                                // ressearching using reverse complement k-mer too
                                reverseComplement(km);
                                direction = 1;
                                find<0, NB_ERR>(delegateParallel, index, km , EditDistance() );
                                
                                // #pragma omp critical
                                // {searchCount[omp_get_thread_num()]++;}
                            }
                        }
                    }
                }
                
                
                keptId.clear();
                // checking number of matches in forward
                for(auto it = results[0].begin(); it != results[0].end(); ++it ){
                    // looking for the longet streak of kmers that are colinear to our read.
                    if(LIS(it->second).size() >= treshold ){
                        keptId[it->first] = it->second;
                        // if read id is associated, do not use it again
                        #pragma omp critical
                        stopSearch.insert(it->first);
                    }
                }

                // exporting
                #pragma omp critical
                printReadNamesWithLIS(r,keptId,ids,outputFile);

                // and reverse if wanted
                if(rc){
                    keptId.clear();
                    
                    for(auto it = results[1].begin(); it != results[1].end(); ++it ){
                        // Reversing the results, otherwise LIS won't really make sens.
                        int_vector reversed = it->second;
                        std::reverse(reversed.begin(), reversed.end());
                        if(LIS(reversed).size() >= treshold ){
                            keptId[it->first] = it->second;
                            // if read id is associated, do not use it again
                            #pragma omp critical
                            stopSearch.insert(it->first);
                        }
                    }
                    #pragma omp critical
                    printReadNamesWithLIS(r,keptId,ids,outputFile);
                }
                keptId.clear();
           }
        }
    }

    if(v>=1){
        print("DONE");
    }
        
}



int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("approxCount");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "input filename"));

    addOption(parser, seqan::ArgParseOption(
        "nt", "nb_thread", "Number of thread to work with",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "v", "verbosity", "Level of details printed out",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "rc", "revComp", "Ressearsh k-mer's rev-comp too "));

    addOption(parser, seqan::ArgParseOption(
        "s", "sampling", "Use read sampling / no reprocess method"));

    addOption(parser, seqan::ArgParseOption(
        "ks", "kmer-skip", "Limit ressearch to 1/ks k-mers in the read",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "nk", "nbKmer", "Minimum number of common k-mer needed to create an edge between two reads.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "k", "k", "Size of the kmers",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "lc", "low-complexity", "Threshold for low complexity",
        seqan::ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, seqan::ArgParseOption(
        "o", "outFile", "path to the output file",
        seqan::ArgParseArgument::STRING, "output file"));

     addOption(parser, seqan::ArgParseOption(
        "i", "index", "path and prefix of the index files",
        seqan::ArgParseArgument::STRING, "index file"));

    // Parse command line.
    
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::string output = "out.edges";     // output file
    std::string indexFile = ""; // index prefix, if any
    unsigned nb_thread = 4;  // default number of thread (4)
    //unsigned nbErr = 2;   // default number of errors EDIT: CAN NOT BE ASSIGNED
    unsigned v = 0;         // verobisty, default = 0
    unsigned k = 30;        // kmerSize, default = 16
    unsigned ks = 3;        // k-mer skip size, default = 3
    unsigned nk= 3;         // Minimum common k-mer, default = 3
    double   lc= 1.25;      // "dust2" low complexity threshold, default = 1.25
    bool rc = isSet(parser, "revComp"); // checking if revcomp is activated
    bool sampling = isSet(parser,"sampling"); // sampling method activation

    getOptionValue(nb_thread, parser, "nt");
    getOptionValue(v, parser, "v");
    getOptionValue(k, parser, "k");
    getOptionValue(nk, parser, "nk");
    getOptionValue(ks, parser, "ks");
    getOptionValue(lc, parser, "lc");
    getOptionValue(output, parser, "o");
    getOptionValue(indexFile, parser, "i");
    
    std::string text;
    getArgumentValue(text, parser, 0);
    
    std::ofstream outputFile;
    outputFile.open (output);
    outputFile << text << "File:";
    outputFile << " nb_thread:" << nb_thread;
    outputFile << " k:" << k;
    outputFile << " nk:" << nk;
    outputFile << " kmer_skipped:" << ks;
    outputFile << " lc:" << lc;
    lc = adjust_threshold( lc, 16, k );
    outputFile << " adusted_lc:" << lc;
    outputFile << " rc:" << rc;
    outputFile << " sampling:" << sampling;
    outputFile << " index:" << indexFile;
    outputFile <<std::endl;
    outputFile <<  " Output:" << output << std::endl;

    approxCount(text,indexFile, k, nb_thread, outputFile,v,ks,nk,rc,lc,sampling);
    outputFile.close();
    if(v>=1){
        print("PROGRAM END");
    }

    return 0;
}