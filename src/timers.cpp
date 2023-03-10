// JEM-Mapper: A C++ implementation for Jaccard Estimate MinHash-based sequence-to-sequence mapping

// Tazin Rahman, Oieswarya Bhowmik, Ananth Kalyanaraman

//      (tazin.rahman@wsu.edu, oieswarya.bhowmik@wsu.edu, ananth@wsu.edu)

// Washington State University

//

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

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

double k_gen_time=0.0, global_k_gen_time=0.0;
double k_shift_time=0.0, global_k_shift_time=0.0;
double kmer_count_time=0.0, global_kmer_count_time=0.0;
double sel_time=0.0, global_sel_time=0.0;
double alltoall_time=0.0, global_alltoall_time=0.0;
double alltoallv_time=0.0, global_alltoallv_time=0.0;
double partial_contig_time=0.0, global_partial_contig_time=0.0;
double preproc_time=0.0, global_preproc_time=0.0;
double gather_time=0.0, global_gather_time=0.0;
double barrier_time=0.0, global_barrier_time=0.0;
double sort_time=0.0, global_sort_time=0.0;
double mod = 0.0, global_mod_time = 0.0;
double sl_time = 0.0, global_sl_time = 0.0;
double hash = 0.0, global_hash_time = 0.0;
double comm_time = 0.0, global_comm_time = 0.0;
double qsl = 0.0, global_qsl_time = 0.0;
double qhash = 0.0, global_qhash_time = 0.0;
double ss = 0.0, global_ss_time = 0.0;
double out = 0.0, global_out_time = 0.0;
double rd = 0.0, global_rd_time = 0.0;
double ag_time = 0.0, global_ag_time = 0.0;
double fl_time = 0.0, global_fl_time = 0.0;

double contigs_time=0.0, global_contigs_time=0.0;
double io_time=0.0, cleanup_time=0.0, vector_time=0.0;
double count_time=0.0, global_count_time=0.0;
double read_input_time=0.0, global_read_input_time=0.0;
double pack_sbuf_time=0.0, global_pack_sbuf_time=0.0;
double unpack_rbuf_time=0.0, global_unpack_rbuf_time=0.0;
double unpack_rbuf_sort=0.0, global_unpack_rbuf_sort=0.0;
double unpack_rbuf_insert=0.0, global_unpack_rbuf_insert=0.0;
double unpack_rbuf_acc=0.0, global_unpack_rbuf_acc=0.0;
double sl_win_time=0.0, global_sl_win_time=0.0;
double sl_win_time_int=0.0, global_sl_win_time_int=0.0;
double sl_lmer_freq=0.0, global_sl_lmer_freq=0.0;

//double kmer_recal_time=0.0, global_kmer_recal_time=0.0;
//double kmer_shift_time=0.0, global_kmer_shift_time=0.0;
double vec_insert_time=0.0, global_vec_insert_time=0.0;
double tmap_insert_time=0.0, global_tmap_insert_time=0.0;


//P2 Timers
double p2alltoall_time=0.0, p2global_alltoall_time=0.0;
double p2alltoallv_time=0.0, p2global_alltoallv_time=0.0;

