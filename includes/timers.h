

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

extern double k_gen_time, global_k_gen_time;
extern double k_shift_time, global_k_shift_time;
extern double kmer_count_time, global_kmer_count_time;
extern double sel_time, global_sel_time;
extern double alltoall_time, global_alltoall_time;
extern double alltoallv_time, global_alltoallv_time;
extern double partial_contig_time, global_partial_contig_time;
extern double preproc_time, global_preproc_time;
extern double gather_time, global_gather_time;
extern double barrier_time, global_barrier_time;
extern double sort_time, global_sort_time;
extern double mod, global_mod_time;
extern double sl_time, global_sl_time;
extern double hash , global_hash_time ;
extern double comm_time , global_comm_time ;
extern double qsl , global_qsl_time ;
extern double qhash , global_qhash_time ;
extern double ss, global_ss_time ;
extern double out, global_out_time;
extern double rd, global_rd_time;
extern double ag_time, global_ag_time;
extern double fl_time, global_fl_time;

extern double contigs_time, global_contigs_time;
extern double io_time, cleanup_time, vector_time;
extern double count_time, global_count_time;
extern double read_input_time, global_read_input_time;
extern double pack_sbuf_time, global_pack_sbuf_time;
extern double unpack_rbuf_time, global_unpack_rbuf_time;
extern double unpack_rbuf_sort, global_unpack_rbuf_sort;
extern double unpack_rbuf_insert, global_unpack_rbuf_insert;
extern double unpack_rbuf_acc, global_unpack_rbuf_acc;
extern double sl_win_time, global_sl_win_time;
extern double sl_win_time_int, global_sl_win_time_int;
extern double sl_lmer_freq, global_sl_lmer_freq;

extern double vec_insert_time, global_vec_insert_time;
extern double tmap_insert_time, global_tmap_insert_time;

//P2 Timers
extern double p2alltoall_time, p2global_alltoall_time;
extern double p2alltoallv_time, p2global_alltoallv_time;

