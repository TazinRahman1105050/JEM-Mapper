# Dependencies:
JEM-Mapper has the following dependencies:

* MPI library (preferably MPI-3 compatible)
* C++14 (or greater) compliant compiler
# Build:
Please specify the k-mer length at the time of building:

For example: make ksize=15


# Execute:
mpiexec -np $number_of_procs $BINARY -s {Contig_Fasta_File} -q {Long_Read_Fasta_File} -a {A_int_Values_File} -b {B_int_Values_File} -p {Prime_int_Values_File} -r $read_segment length -n $NO_OF_TRIALS

For example:

mpiexec -np 4 ./jem -s ~/Ecoli_reads_100x_contigs.fasta -q ~/Ecoli_reads_10x_long_reads.fasta -a ~/A.txt -b ~/B.txt -p ~/Prime.txt -r 1000 -n 30

Notes:
* We have run our codes in a high-performance computing cluster (HPC) Aeolus. FOr this system we had to set the number of processes in the given way. If you use any other system please change the parameters accordingly

# Input arguments meaning
1. -s: input contigs fasta file
2. -q: input long reads fasta file
3. -a: For diffent trials we have used a linear congruential hash funtion of the form: (Ax+B)%P, whwere A and B are integers and P is a prime and x is a kmer we want to hash. SO, using this parameter, we provide the input values for different A values
4. -b: input B values
5. -p: input prime numbers
6. -r: read segment length
7. -n: number of trials
