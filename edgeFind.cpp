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
using index_t = Index<DnaString, BidirectionalIndex<FMIndex<void,TFastConfig> > >;


const auto boot_time = std::chrono::steady_clock::now();
// Shortcut to print text in stdout
template<typename TPrintType>
void print(TPrintType text)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" << text << std::endl;
}


inline  bool haveLowComplexity(DnaString sequence, double threshold){
    // New version, using "DUST2" method
    // scanning 2-mers, squaring count, discard if over limit
    
    unsigned l = length(sequence);
    std::unordered_map<std::string,int> counter;
    std::string seq = toCString(CharString(sequence));
    // reading using sliding window of 2
    for(int i =0; i < l-1; i++){
        std::string c = seq.substr(i,2);
        if(counter.count(c) !=0 ){
            counter[c]+=1;    
        }
        else{
            counter[c]=1;
        }
    }
    float s = 0;
    float sum = 0.0;
    for(auto v:counter){
        sum+=  float(v.second * (v.second-1) / 2);  
    }
    s =  sum / (l-2);
    //std::cout << "LOW C "<<sequence << " "  << s << " "<< threshold << " " << bool(s>=threshold) << std::endl;
    if(s>= threshold){
        
        return(true);
    }
    return(false);
}



template<typename TNumberType, typename TArrayType> 
inline int dichoFind(TNumberType target, TArrayType &values){
    int up = values.size();
    
    int down = 0;
    int mid;

    while(up-1!=down){
        mid = (up+down)/2;
        if(values[mid]>target){
            up = mid;
        }
        else if(values[mid]<=target){
            down = mid;
        }

    }   
    return(down);
}


// Function used to print the "clustering" results
void printReadNames(int currentReadId, read2pos_map_t readIdMap , StringSet<CharString> realIds ,std::ofstream &outputFile)
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

    outputFile << std::endl;
}


// Search and count a kmer list in a fasta file, at at most a levenstein distance of 2.
void approxCount(const std::string& filename, const std::string & indexFile, const int k, const int& nbThread, std::ofstream &outputFile, int v, int bloc_size, int treshold, bool rc, double lc, bool sampling, int ob){
    
    if(v>=1)
        print("PARSING FILE");

    // Parsing input fasta file
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn(toCString(filename));
    readRecords(ids, seqs, seqFileIn);

    // Creating the sequence to index
    DnaString allReads;
    std::vector<unsigned long int> lenReads; // summed size of all the reads
    int n = 0;
    lenReads.push_back(n);
    for (unsigned i = 0; i < length(ids); ++i){
        n += length(seqs[i]);
        lenReads.push_back(n);
        append(allReads,seqs[i]);
    }
    if(v>=1)
        print("DONE");

    // Constants (should be in upper case)
    const int nbRead  = length(ids);     // Number of reads in the file

    index_t index;

    // INDEX CREATION ----------
    if(v>=1)
        print("CREATING INDEX");
        if(open(index,indexFile.c_str())){
            print("LOADED FROM INDEX FILE");
        }
        else{
            print("NO INDEX FILE, CREATING FROM SCRATCH...");
            index_t index(allReads);
        }
    if(v>=1)
        print("DONE");
    // -------------------------

    // RESSEARCH ---------------
    omp_set_num_threads(nbThread);
    if(v>=1)
        print("STARTING RESSEARCH");
    std::set<read_id_t> stopSearch;
    int searchCount = 0;
    #pragma omp parallel  firstprivate(index) shared(ids, seqs, treshold, rc, stopSearch, lc, sampling, searchCount)
    {
        // storing results
        std::array<read2pos_map_t, 2> results;
        int direction = 0 ;
        read2pos_map_t keptId;
        std::unordered_map<std::string, bool> processedKmer;
        StringSet<DnaString> kmSet;

        // Max number of errors
        const uint8_t NB_ERR = 1;

        auto delegateParallel = [&](auto & iter, const DnaString & needle, int errors)// __attribute__((noinline))
        {
            
            for (auto occ : getOccurrences(iter)){
                
                // Identifying read Id and position on read    
                // and esting if the occurence is on one read only
                // (not between two reads)
                int readId= dichoFind(occ,lenReads);
                int l = lenReads[readId+1] - lenReads[readId];
                read_pos_t readPos =  occ - lenReads[readId];
                                
                if( readPos + k < l )
                {   
                    #pragma omp critical
                    if(std::find(results[direction][readId].begin(), results[direction][readId].end(), readPos/ob ) == results[direction][readId].end() ){
                        results[direction][readId].push_back(readPos/ob);
                    }
                }
            }
        };

        bool not_in;
        // going through the read list
        
        #pragma omp for schedule(dynamic)
        for(read_id_t r=0; r<nbRead; r++)
        //for(int r=0; r<1; r++)
        {   
            if(v>=1 and (r+1) %(nbRead/100) == 0){
                print(float(r+1)/nbRead*100);

            }

            // checking read id is not already found, if sampling is activated
            if(sampling){
                not_in = stopSearch.find(r) == stopSearch.end();
            }
            else{
                not_in = true;
            }
            // if ressearsh is allowed, then proceed
            if(not_in){

                // cleaning result sets and temporary structures
                results[0].clear();
                results[1].clear();
                processedKmer.clear();
                keptId.clear();
                clear(kmSet);
                DnaString readSequence = seqs[r];
                
                // going through the read using K size window
                for( int i = 0; i < length(readSequence) - k + 1; i+= bloc_size ){
                //for( int i = 0; i < 1; i++ ){
                    Infix<DnaString>::Type inf = infix(readSequence, i, i+k);
                    DnaString km = DnaString(inf);
                    std::string strKm = toCString(CharString(km));
                    //print(strKm);
                    if(not processedKmer[strKm]){
                        processedKmer[strKm] = true;
                        if(not haveLowComplexity(km,lc)){
                            appendValue(kmSet, km);
                        }
                    }
                }

                direction = 0;
                find<0, NB_ERR>(delegateParallel, index, kmSet , EditDistance(), Parallel());
                searchCount++;
                // if needed:
                if(rc){
                    for(auto& el: kmSet){
                        reverseComplement(el);
                    }
                    
                    direction = 1;
                    find<0, NB_ERR>(delegateParallel, index, kmSet , EditDistance(), Parallel() );
                    
                }
                
                // checking number of matches in forward
                for(auto it = results[0].begin(); it != results[0].end(); ++it ){
                    if(it->second.size() >= treshold ){
                        keptId[it->first] = it->second;
                        // if read id is associated, do not use it again
                        if(sampling){
                            #pragma omp critical
                            stopSearch.insert(it->first);
                            searchCount++;
                        }
                    }
                }

                // exporting
                printReadNames(r,keptId,ids,outputFile);

                // and reverse if wanted
                if(rc){
                    keptId.clear();
                    for(auto it = results[1].begin(); it != results[1].end(); ++it ){
                        if(it->second.size() >= treshold ){
                            keptId[it->first] = it->second;
                            // if read id is associated, do not use it again
                    
                            if(sampling){
                               #pragma omp critical
                               stopSearch.insert(it->first);

                            }
                        }
                    }
                    
                    printReadNames(r,keptId,ids,outputFile);
                }


            #pragma omp critical
            if(rc){
                searchCount++;
            }
            #pragma omp critical
            searchCount++;
            }
           
        }
        
    }
    std::cout << searchCount << " ressearchs have been made. (RC included)" << std::endl;
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
        "nt", "nbThread", "Number of thread to work with",
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
        "ob", "occurence-bloc", "compact occurences in blocs of size ob (11,50,75 are in the same bloc for ob=100",
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
    std::string indexFile = "no_index"; // index prefix, if any
    unsigned nbThread = 4;  // default number of thread (4)
    //unsigned nbErr = 2;     // default number of errors EDIT: CAN NOT BE ASSIGNED
    unsigned v = 0;         // verobisty, default = 0
    unsigned k = 30;        // kmerSize, default = 16
    unsigned ks = 3;        // k-mer skip size, default = 3
    unsigned nk= 3;         // Minimum common k-mer, default = 3
    double lc= 1.25;        // "dust2" low complexity threshold, default = 1.25
    unsigned ob= 100;       // occurence bloc default fold size = 100
    bool rc = isSet(parser, "revComp"); // checking if revcomp is activated
    bool sampling = isSet(parser,"sampling"); // sampling method activation

    getOptionValue(nbThread, parser, "nt");
    getOptionValue(v, parser, "v");
    getOptionValue(k, parser, "k");
    getOptionValue(nk, parser, "nk");
    getOptionValue(ks, parser, "ks");
    getOptionValue(ob, parser, "ob");
    getOptionValue(lc, parser, "lc");
    getOptionValue(output, parser, "o");
    getOptionValue(indexFile, parser, "i");
    
    std::string text;
    getArgumentValue(text, parser, 0);
    
    std::ofstream outputFile;
    outputFile.open (output);
    if(v>=1){
        outputFile << text << " File:" << text << " nbThread:" << nbThread << " k:" << k;
        outputFile <<  " nk:" << nk << " ob:" << ob << " lc:" << lc;
        outputFile <<  " rc:" << rc << " sampling:" << sampling <<std::endl;
        outputFile <<  " Output:" << output << std::endl;
    }

    approxCount(text,indexFile, k, nbThread, outputFile,v,ks,nk,rc,lc,sampling,ob);
    outputFile.close();
    if(v>=1){
        print("END");
    }

    return 0;
}