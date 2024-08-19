We present an efficient parallel algorithmic workflow, called JEM-mapper, that uses a new minimizer-based Jaccard estimator (or JEM) sketch to perform alignment-free mapping of long reads. We are able to implement and test this new multithreaded extension to the code (MPI+OpenMP multithreading). 
![image](https://github.com/user-attachments/assets/e0df386c-64a4-4237-800d-648af513b4e1)


# Citation information:
Rahman, Tazin, Oieswarya Bhowmik, and Ananth Kalyanaraman. "An Efficient Parallel Sketch-based Algorithm for Mapping Long Reads to Contigs." 2023 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW). IEEE, 2023. DOI: 10.1109/IPDPSW59300.2023.00037

Rahman, Tazin, Oieswarya Bhowmik, and Ananth Kalyanaraman. "An Efficient Parallel Sketch-based Algorithmic Workflow for Mapping Long Reads." bioRxiv (2023): 2023-11.

# Dependencies:
JEM-Mapper has the following dependencies:

* MPI library (preferably MPI-3 compatible)
* C++14 (or greater) compliant compiler

# Build:
make ksize=$KMER_SIZE

For example:
make ksize = 15

# Execute:
export OMP_NUM_THREADS= $number_of_threads     
mpiexec -np $number_of_procs $BINARY -s {Contig_Fasta_File} -q {Long_Read_Fasta_File} -a {A_int_Values_File} -b {B_int_Values_File} -p {Prime_int_Values_File} -r $read_segment length -T $NO_OF_TRIALS

For example: if we want to run it on 8 threads and 4 processes:  
export OMP_NUM_THREADS=8  
mpiexec -np 4 ./jem -s ~/Ecoli_reads_100x_contigs.fasta -q ~/Ecoli_reads_10x_long_reads.fasta -a ~/A.txt -b ~/B.txt -p ~/Prime.txt -r 1000 -T 30

Notes:
* This code has been tested on high-performance computing cluster (HPC) with MPI compatibility. For the system we used we had to set the number of processes in the given way. Please change the parameters accordingly.

# Input arguments 
* -s: input contigs fasta file
* -q: input long reads fasta file
* -a: For diffent trials we have used a linear congruential hash funtion of the form: (Ax+B)%P, whwere A and B are integers and P is a prime and x is a kmer we want to hash. So, using this parameter, we provide the input values for different A values
* -b: input B values
* -p: input prime numbers
* -r: read segment length
* -n: number of trials
