#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"
#include "timers.h"
#include <iostream>
#include <fstream>
#include <sstream>
//#include "serial.h"

using std::cout; using std::cerr;
using std::endl; using std::string;
using std::ifstream; using std::ostringstream;

long int MAX_KMER_COUNT=0;
int rank, size;
int coverage=0;
std::string inputFileName;
std::string queryFileName;
int read_length=0;
int num_buckets=0;
int num_threads=0;
int num_contigs=0;
int node_threashold=0;

int num_batch_transfers=0;
long int this_contig_id = 0;

std::vector<lmer_t> lmer_frequency(LMER_SIZE,0);
std::vector<lmer_t> global_lmer_frequency(LMER_SIZE,0);

/* Vector containing k-mers allocated to each process */
std::vector<KmerPairs> kmer_proc_buf;
std::vector<std::vector<kmer_t>> kmer_sets;
//std::vector<std::unordered_map<kmer_t, std::vector<int>> > Tl;

//kmer_t *hasha, *hashb;
uint64_t hasha = 68111;
uint64_t hashb = 105929;




void parseCommandLine(const int argc, char * const argv[]);

int retrieve_proc_id (lmer_t min_lmer)
{
    //determine which proc holds the l-mer bucket based on a block-cyclic
    //distribution of buckets across the processes

    //int proc_num = (int)((int)(floor(min_lmer/4)) % (int)size); // static block-cyclic
    //int proc_num = (int)((int)(hash31(hasha[min_lmer], hashb[min_lmer], min_lmer)) % (int)size); // randomized
    int proc_num = (int)((int)(hash31(hasha, hashb, min_lmer)) % (int)size); // randomized

    return proc_num;
}
#ifdef EXTEND_KMER
int retrieve_proc_id (kmer_t min_lmer)
{
    //determine which proc holds the l-mer bucket based on a block-cyclic
    //distribution of buckets across the processes

    //int proc_num = (int)((int)(floor(min_lmer/4)) % (int)size); // static block-cyclic
    //int proc_num = (int)((int)(hash31(hasha[min_lmer], hashb[min_lmer], min_lmer)) % (int)size); // randomized

    //int proc_num = (int)((int)(hash31(hasha, hashb, min_lmer)) % (int)size); // randomized
    int proc_num = (int)((int)(uhash31(hasha, hashb, min_lmer)) % (int)size);

    return proc_num;
}
#endif

void set_num_threads()
{

#pragma omp parallel
{
    num_threads = omp_get_num_threads();
}

}

std::string readFileIntoString(const std::string& path) {
    ifstream input_file(path);
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '"
             << path << "'" << endl;
        exit(EXIT_FAILURE);
    }
    return string((std::istreambuf_iterator<char>(input_file)), std::istreambuf_iterator<char>());
}

int get_file_size(std::string filename) // path to file
{
    FILE *p_file = NULL;
    p_file = fopen(filename.c_str(),"rb");
    fseek(p_file,0,SEEK_END);
    int size = ftell(p_file);
    fclose(p_file);
    return size;
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    parseCommandLine(argc, argv);
    set_num_threads();

    //std::string contigFileName = inputFileName + std::to_string(rank) + ".fa";
    //std::cout<<std::to_string(rank)<<" - "<<inputName<<"\n";
    double time_l1 = MPI_Wtime ();
    
    //input_read_data cdata = perform_input_reading(rank, size, inputFileName, 108539);
   // input_read_data cdata = perform_input_reading(rank, size, inputFileName, 25980); //34680 //25980 EC //27586 PA
    input_read_data cdata = perform_input_reading(rank, size, inputFileName, 127978);
    //input_read_data rdata = perform_input_readingC(rank, size, inputName, read_length);
    
    //string file_contents;
    //size_t file_size;
    
    //std::cout<<"Print statement\n";
    //std::cout<<contigFileName<<"\n";
        //contigdata = perform_input_readingC(rank, size, contigFileName, read_length);
    //file_contents = readFileIntoString(contigFileName);
    //file_size = get_file_size(contigFileName);
    //std::cout<<"Print statement\n";
        
        //file_contents = readFileIntoString(filename);
    
    //read_array();
    //perform_kmer_counting (file_contents, file_size);
    int M;
    int total_subjects;
    
    
    input_read_data rdata = perform_input_reading(rank, size, queryFileName, read_length);
    double time_l3 = MPI_Wtime ();
    generate_set_of_subjects (cdata.read_data, cdata.read_data_size, cdata.start_index, rdata.read_data, rdata.read_data_size, rdata.start_index, &M, &total_subjects);
    //int total_number_of_subs_in_p = kmer_sets.size();
    //kmer_t** Hash_table = new kmer_t*[total_subjects];
    
    /*int total_hash_functions = 150;
    for (int i = 0; i < total_subjects; i++) {
        Hash_table[i] = new kmer_t[total_hash_functions];
    }
    genereate_hash_table(M, total_subjects, Hash_table);*/
    //double time_l11 = MPI_Wtime ();
    
    double time_l2 = MPI_Wtime ();
    double f_time = time_l2 - time_l3;
    //printf ("%d Average time for out across all procs (secs): %f \n", rank, f_time);
    //free(cdata.read_data);
    //free(rdata.read_data);
    /*
    if(rank == 0)
    {
        const std::string s0("/home/trahman/GitRepo/pakman-pakman_latest/Hash_");
        char proc_id[3];
        char output_file_name[25];
        sprintf(proc_id, "%d", rank); 
    //strcpy(output_file_name,"wired_mn_out_");
        strcpy(output_file_name, s0.c_str());
        strcpy(&output_file_name[strlen(output_file_name)],proc_id);
        strcpy(&output_file_name[strlen(output_file_name)],".log");
        FILE *f = fopen(output_file_name, "w");
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        for(int h = 0; h<total_subjects; h++)
        {
            for(int k = 0; k<200; k++)
            {
                //fprintf(f,"%ld ", Hash_table[h][k]); 
            }
        //    fprintf(f, "\n");
        }
    }
    
    */
    //double time_l3 = MPI_Wtime ();
    //input_read_data rdata = perform_input_reading(rank, size, queryFileName, read_length);


    //perform_kmer_counting (rdata.read_data, rdata.read_data_size);
   // generate_set_of_queries (rdata.read_data, rdata.read_data_size, rdata.start_index, total_subjects, M, Hash_table);
    
    double time_l4 = MPI_Wtime ();
    double s_time = time_l4 - time_l3;
    //printf ("%d Average time for l across all procs (secs): %f \n", rank, s_time);
#ifndef PERFORM_CONTIG_GENERATION

    free_kmer_count_buffers();

#endif

#ifdef PERFORM_CONTIG_GENERATION

    /************ CONTIG GENERATION BEGINS *******************/
    /* Phase 1) Macro Node construction and Initial Wiring
     * Phase 2) Construction of Independent Set. Merge and delete all candidate nodes from Independent set.
     * Phase 3) Gather all remaining nodes
     * Phase 4) Enumerate/Generate the contigs with the Walk algorithm
     */


    /* PHASE 1: MACRO NODE CONSTRUCTION */
//    std::vector<std::pair<kmer_t,MacroNode>> MN_map;
//    begin_mnode_construction(MN_map);


    /* PHASE 1: MACRO NODE WIRING */
//    initiate_mnode_wiring(MN_map);

#ifdef DEBUG_WIRE_INIT
    debug_wired_mnodes(MN_map);
#endif

    /* PHASE 2: INDEPENDENT SET CONSTRUCTION */
    
    // temp vector for storing all partial contigs generated during Phase 2
//    std::vector<BasePairVector> partial_contig_list;
//    size_t global_num_nodes = begin_iterative_compaction(MN_map, partial_contig_list);

//    MPI_Barrier(MPI_COMM_WORLD);

#ifdef COMPACT_PGRAGH
    //retain a list of terminal prefixes for each individual process, potential begin k-mers
  //  std::vector<BeginMN> list_of_begin_kmers;
    //identify_begin_kmers (MN_map, list_of_begin_kmers);

#ifdef DEBUG_BKMERS
    //debug_begin_kmers_list (MN_map, list_of_begin_kmers);
#endif

    /* Perform an Allgather such that all macro_nodes are accessible to all procs */

    //std::vector<std::pair<kmer_t,MacroNode>> global_MN_map(global_num_nodes); // new map for storing all the macro_nodes
    //generate_compacted_pakgraph(MN_map, global_MN_map);

    //traverse_pakgraph(global_MN_map, list_of_begin_kmers, partial_contig_list);
#endif // end of COMPACT_PGRAGH

#endif // end of PERFORM_CONTIG_GENERATION

   // MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

void parseCommandLine(const int argc, char * const argv[])
{
  int ret;

  while ((ret = getopt(argc, argv, "s:q:b:r:c:t:n:")) != -1) {
    switch (ret) {
    case 's':
       inputFileName.assign(optarg);
       //std::cout << inputFileName << std::endl;
       break;
    case 'q':
       queryFileName.assign(optarg);
       //std::cout << inputFileName << std::endl;
       break;
    case 'b':
       MAX_KMER_COUNT = atol(optarg);
       //std::cout << MAX_KMER_COUNT << std::endl;
       break;
    case 'r':
       read_length = atoi(optarg);
       //std::cout << read_length << std::endl;
       break;
    case 'c':
       coverage = atoi(optarg);
       //std::cout << coverage << std::endl;
       break;
    case 't':
       num_buckets = atoi(optarg);
       //std::cout << num_buckets << std::endl;
       break;
    case 'n':
       node_threashold = atoi(optarg);
       //std::cout << node_threashold << std::endl;
       break;
    default:
       assert(0 && "Should not reach here!!");
       break;
    }
  }

  if (rank==0 && inputFileName.empty()) {
      std::cerr << "Must specify an input FASTA file name with -f" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

//  if (rank==0 && !MAX_KMER_COUNT) {
//      std::cerr << "Must specify a batch size for k-mer counting with -b" << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, -99);
//  }

  if (!MAX_KMER_COUNT) {
      MAX_KMER_COUNT=100000000;
      if (rank==0) std::cout << "batch size not specified with -b, set to default value MAX_KMER_COUNT=100000000" << std::endl;
  }

  if (rank==0 && !read_length) {
      std::cerr << "Must specify read_length with -r" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

  if (rank==0 && !coverage) {
      std::cerr << "Must specify coverage of the input data with -c" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
/*
  if (rank==0 && (read_length>250 || read_length<100)) {
      std::cerr << "Must provide short reads of length >100 and <=250 with -r" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }*/

  if (!num_buckets) {
      num_buckets=21;
      if (rank==0) std::cout << "num_buckets not specified with -t, set to default value num_buckets=21" << std::endl;
  }

  if (!node_threashold) {
      node_threashold=100000;
      if (rank ==0 )std::cout << "node_threashold not specified with -n, set to default value node_threashold=100K" << std::endl;
  }

  if (rank == 0) {
           printf("K-mer size: %d, L-mer size: %d, Number of Processes: %d, MAX_KMER_COUNT: %ld Coverage: %d,",
           WINDW_SIZE+1, LMER_LENGTH, size, MAX_KMER_COUNT, coverage);
           printf("Avg read length: %d, Num buckets: %d, Node threshold: %d\n", 
           read_length, num_buckets, node_threashold);
  }


} // parseCommandLine
