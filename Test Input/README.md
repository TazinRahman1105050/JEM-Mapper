# Generate contigs using Minia assembler:

To generates contigs from short reads we have used minia short reads assembler (https://github.com/GATB/minia).
Please find their complete manual: https://github.com/GATB/minia/raw/master/doc/manual.pdf

The basic usage is as follows where input_file is the input short read fasta file:
./minia -in [input file] -kmer-size [kmer size] -out [prefix]


Input long read files are generated using SimIt simulator
