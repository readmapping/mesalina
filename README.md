MESALINA
==========

Repository for the program mesalina.

mesalina
----------

To build mesalina, just execute `./make` .

To run mesalina, we added some example data.
First we need to build the index:

	./mesalina index -s <sparseness> -x <reference.fasta> -p <index name>

For example:

	./mesalina index -s 1 -x demo/TAIR10_chr1.fas -p demo/TAIR10_chr1

Now we can align reads to the reference:

	./mesalina splice -S <output> -i <index name> -x <reference.fasta> -U <queries.fasta>

For example:

	./mesalina splice -S demo/output.sam -i demo/TAIR10_chr1 -x demo/TAIR10_chr1.fas -U demo/AT1G01530.1.fa
	
When you want to view the alignment, add the parameter `-A`:
	
	./mesalina splice -A -S demo/output.sam -i demo/TAIR10_1_sparse-1 -x demo/TAIR10_chr1.fas -U demo/AT1G01530.1.fa
	