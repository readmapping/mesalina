GSLAP
==========

Repository for the program GSLAP.

GSLAP
----------

To build gslap, just execute `./make` .

To run gslap, we added some example data.
First we need to build the index:

	./gslap index -s <sparseness> -x <reference.fasta> -p <index name>

For example:

	./gslap index -s 1 -x demo/TAIR10_chr1.fas -p demo/TAIR10_chr1

Now we can align reads to the reference:

	./gslap splice -S <output> -i <index name> -x <reference.fasta> -U <queries.fasta>

For example:

	./gslap splice -S demo/output.sam -i demo/TAIR10_chr1 -x demo/TAIR10_chr1.fas -U demo/AT1G01530.1.fa
	
When you want to view the alignment, add the parameter `-A`:
	
	./gslap splice -A -S demo/output.sam -i demo/TAIR10_1_sparse-1 -x demo/TAIR10_chr1.fas -U demo/AT1G01530.1.fa
	