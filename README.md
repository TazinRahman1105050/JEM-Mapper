# Dependencies:
JEM-Mapper has the following dependencies:

* MPI library (preferably MPI-3 compatible)
* C++14 (or greater) compliant compiler
# Build:
Please specify the k-mer length at the time of building:

For e.g: make ksize=15


# Execute:
mpiexec -np $procs $BINARY -s $INPUT_SUBJET -q $INPUT_QUERY -n $NO_OF_TRIALS

For example:

mpiexec -np 4 ./jem -s ~/Ecoli_reads_100x_contigs.fasta -q ~/Ecoli_reads_10x_long_reads.fasta -n 30
