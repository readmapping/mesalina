GSLAP
==========

Repository for the program GSLAP.
The program name is currently changing from Al-fisfisa (A Long Fragment-Intron Segment-Fragment-Intron Segment Aligner) to GSLAP (Genomic Spliced Long-read Alignment Program). All references to the name Al-fisfisa will soon be replaced by GSLAP.

Al-fisfisa (older name for GSLAP)
----------

To build alfisfisa, just execute `./make` .

To run alfisfisa, we added some example data.
First we need to build the index:

	./alfisfisa index -s <sparseness> -x <reference.fasta> -p <index name>

For example:

	./alfisfisa index -s 1 -x demo/TAIR10_chr1.fas -p demo/TAIR10_chr1

Now we can align reads to the reference:

	./alfisfisa splice -S <output> -i <index name> -x <reference.fasta> -U <queries.fasta>

For example:

	./alfisfisa splice -S demo/output.sam -i demo/TAIR10_chr1 -x demo/TAIR10_chr1.fas -U demo/AT1G01530.1.fa
	
When you want to view the alignment, add the parameter `-A`:
	
	./alfisfisa splice -A -S demo/output.sam -i demo/TAIR10_1_sparse-1 -x demo/TAIR10_chr1.fas -U demo/AT1G01530.1.fa
	