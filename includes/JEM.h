
// JEM-Mapper: A C++ implementation for Jaccard Estimate MinHash-based sequence-to-sequence mapping

// Tazin Rahman, Oieswarya Bhowmik, Ananth Kalyanaraman

//      (tazin.rahman@wsu.edu, oieswarya.bhowmik@wsu.edu, ananth@wsu.edu)

// Washington State University

//
//For citation, please cite the following paper:
//An efficient parallel sketch-based algorithm for mapping long reads to contigs.
//Tazin Rahman, Oieswarya Bhowmik, Ananth Kalyanaraman.
//Proc. 2023 IEEE International Workshop on High Performance Computational Biology (HiCOMB)

// **************************************************************************************************

// Copyright (c) 2023. Washington State University ("WSU"). All Rights Reserved.
// Permission to use, copy, modify, and distribute this software and its documentation
// for educational, research, and not-for-profit purposes, without fee, is hereby
// granted, provided that the above copyright notice, this paragraph and the following
// two paragraphs appear in all copies, modifications, and distributions. For
// commercial licensing opportunities, please contact The Office of Commercialization,
// WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
// commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

// IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
// OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
// THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
// ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
// OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
// **************************************************************************************************

#ifndef DISTRIBUTE_KMERS_H
#define DISTRIBUTE_KMERS_H

#include <stdio.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <queue>          // std::priority_queue
#include <string>
#include <assert.h>
#include <cstring>
#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <atomic>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include<cstdint>


#define KMER_LENGTH                    (WINDW_SIZE+1)
#define LMER_SIZE                      (pow(4, LMER_LENGTH))
#define MN_LENGTH                      (KMER_LENGTH-1)

#define MOD 2147483647
#define HL 31


#ifdef EXTEND_KMER
//static_assert(std::is_same_v<__uint128_t, unsigned __int128>);
//typedef unsigned __int128 uint128_t;
//typedef uint128_t kmer_t;
typedef __uint128_t kmer_t;
#else
typedef uint64_t kmer_t;
#endif

using BasePair = uint8_t;
typedef BasePair ElType;
typedef uint64_t lmer_t;



#define KMER_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*KMER_LENGTH))
#define BASE_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*32))
#define LMER_MASK           ((~0UL)) >> ((sizeof(lmer_t)*8) - (2*16))
#define MN_MASK             ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*MN_LENGTH))
#define SUCC_MASK           ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*(KMER_LENGTH-1)))








/*struct KeyEqual
{
      std::size_t operator()(const BasePairVector& k1, const BasePairVector& k2) const
      {
          kmer_t kmer1 = 0, kmer2 = 0;
          for (kmer_t i=0; i<k1.size(); i++) {
               kmer1 = mnmer_shift(kmer1, k1[i]);
          }
          for (kmer_t i=0; i<k2.size(); i++) {
               kmer2 = mnmer_shift(kmer2, k2[i]);
          }
          
          if (kmer1 == kmer2)
              return true;
          else
              return false;
      }
};
*/

class Comp_rev{
    const std::vector<std::pair<int, int>> & _v;
  public:
    Comp_rev(const std::vector<std::pair<int, int>> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return ((_v[i].second > _v[j].second) ||
                 ((_v[i].second == _v[j].second) &&
                  ( _v[i].first > _v[j].first))
                );

   }
};

class Comp_rev1{
    const std::vector<std::pair<kmer_t, int>> & _v;
  public:
    Comp_rev1(const std::vector<std::pair<kmer_t, int>> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return ((_v[i].first < _v[j].first) ||
                 ((_v[i].first == _v[j].first) &&
                  ( _v[i].second > _v[j].second))
                );

   }
};




typedef struct __attribute__ ((__packed__)) MinHash_pairs
{ 
  int trail;
  kmer_t seq;
  int subject_id;

} MinHashPairs;
static_assert(sizeof(MinHashPairs) == sizeof(int)+(sizeof(kmer_t)+sizeof(int)), "struct MinHashPairs shouldn't be padded");

typedef struct __attribute__ ((__packed__)) Top_Hit
{ 
  int sub;
  int score;

} TopHit;
static_assert(sizeof(MinHashPairs) == sizeof(int)+(sizeof(kmer_t)+sizeof(int)), "struct MinHashPairs shouldn't be padded");

//data structure for storing Reads
typedef struct Rd_Sequence
{
	// Read data
	char *read_data;
	size_t read_data_size;
  int start_index;

} input_read_data;




class Comp{
    const std::vector<kmer_t> & _v;
  public:
    Comp(const std::vector<kmer_t> & v) : _v(v) {}
    bool operator()(size_t i, size_t j){
         return _v[i] < _v[j];
   }
};

inline ElType kmerel(kmer_t kmer, unsigned pos) {
  assert(pos < KMER_LENGTH);

  return ElType(((kmer) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3));
}

inline ElType lmerel(lmer_t kmer, unsigned pos) {
  assert(pos < LMER_LENGTH);

  return ElType(((kmer) >> ((LMER_LENGTH-(pos)-1)*2)) & (0x3));
}

inline ElType kmerel_mn(kmer_t kmer, unsigned pos) {
      assert(pos < MN_LENGTH);

        return ElType(((kmer) >> ((MN_LENGTH-(pos)-1)*2)) & (0x3));
}


inline lmer_t kmer_to_lmer(kmer_t kmer_in, unsigned pos, lmer_t kmer_out) {
  assert(pos < KMER_LENGTH);

  //ElType int_el = ElType(((kmer_in) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3));
  //assert(int_el>=A && int_el<=G);
  //return (kmer_t)((kmer_out<<2) | (kmer_t)int_el) & (kmer_t)LMER_MASK;

  return (lmer_t)((kmer_out<<2) | (lmer_t)(ElType(((kmer_in) >> ((KMER_LENGTH-(pos)-1)*2)) & (0x3)))) & (lmer_t)LMER_MASK;
}

inline ElType char_to_el(char ch) {
              return (ElType)((((ElType)ch)>>1) & 0x7);
}


inline kmer_t kmer_cons(kmer_t kmer_in, 
                        unsigned pos,
                        ElType el) {
  //assert(el>=A && el<=G);
  assert(pos < KMER_LENGTH);
 
  return kmer_in | ((kmer_t)el << ((KMER_LENGTH-pos-1)*2));
}

inline kmer_t kmer_shift(kmer_t kmer_in, 
                        ElType el) {

  //assert(el>=A && el<=G);
  
  return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)KMER_MASK;
  //return ((kmer_in<<2) | (kmer_t)el) & KMER_MASK;
}

inline lmer_t lmer_shift(lmer_t kmer_in, 
                        ElType el) {

  //assert(el>=A && el<=G);
  
  return (lmer_t)((kmer_in<<2) | (lmer_t)el) & (lmer_t)LMER_MASK;
  //return ((kmer_in<<2) | (kmer_t)el) & KMER_MASK;
}




 

void parse_alphabet_file (FILE * fp);
void free_reads (int num_reads);

void change_to_num (char *pred_kmer, int *num);
void change_to_char (char *pred_kmer, int num);
int find_max (int array[]);
int find_min (int array[], int *min);
int compare (const void * a, const void * b);
int FindIndex( const int a[], int size, int value);
char* DivideReads(MPI_File *in, const int rank, const int size, 
        const int overlap, uint64_t *nlines, size_t *data_size);
//input_read_data perform_input_reading (const int rank, 
//        const int size, char *argv[]);
input_read_data perform_input_reading (const int rank, const int size,
                                       std::string &fileName, int read_length);
void perform_kmer_counting (std::string read_data1, size_t length1);
void print_kmers (int klen, int min_t);
//void SortAndAggregate(std::vector<kmer_t>& arr, std::vector<int>& count);

void Sliding_window_l (const char *ptr, size_t length);
void Sliding_window (char *ptr, size_t length, int *M_for_individual_process, int *num_subjects,
                     std::vector<MinHashPairs> &initial_sets, int s_index);
//void Sliding_window (const char *ptr, size_t length, int *n_kmers, 
//                     std::vector<std::vector<kmer_t>> &kmers_per_proc, 
//                     std::vector<std::vector<int>> &kmer_cnt_tmp_buf);
//void transfer_kmers (std::vector<int>& scounts_kmer, 
//                     std::vector<std::vector<KmerPairs>> &kmers_per_proc); 








//std::vector<ModNodeInfo> serialize_and_transfer 
//                         (std::vector< std::vector<ModNodeInfo> >& mn_nodes_per_proc);
/*void serialize_and_transfer 
                         (std::vector< std::vector<ModNodeInfo> >& mn_nodes_per_proc,
                          std::vector<std::pair<kmer_t,MacroNode>> &MN_map,
                          std::vector<size_t> &rewire_pos_list,
                          int num_itr);
*/



void generate_modified_set(int M, std::vector<std::vector<kmer_t>> &previous_sets, int N);

int convert_to_int (char given_char);
input_read_data perform_input_readingC (const int rank, const int size, std::string &fileName, int read_length);
std::string readFileIntoString(const std::string& path);
int get_file_size(std::string filename);
void read_array();
void get_hash_value(kmer_t **A1, int M, kmer_t **Prefix);
void generate_set_of_subjects (char *ptr, size_t length, int s_index,char *read_data, size_t r_length, int start_index, int *M_final, int *num_subjects);
void genereate_hash_table(int M, int total_subjects, kmer_t **Ag_Hash_Table);
void get_hash_value_queires(std::vector<std::vector<kmer_t>> &modified_sets, kmer_t **A1);
void Sliding_window_queires (char *ptr, size_t length, int *num_queries,
                     std::vector<std::unordered_map<kmer_t, std::vector<int>> > Tl, int start, int total_subjects);
void generate_set_of_queries (const char *read_data, size_t length, int start_index, int total_subjects, int M, kmer_t **Ag_Hash_Table);
void generate_modified_set_queries(int M, std::vector<std::vector<kmer_t>> &previous_sets);
char convert_to_char (char given_char);

//////////////////////////////////////////////////////////////
#endif
