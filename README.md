# Origin-of-Replication-Finder
Bioinformatic Algorithm to Find the Origin of Replication on Bacterial Genomes

It relies on the fact on the biomechanics of bacterial division, where GC-skew reaches a minimum around the origin of replication. It scans for consensus sequences around this area of the genome.

This can be run from the terminal like so as 'python main.py -f Genome_E_coli.fas -m 1 -k 9 -w 1000'. An example E coli genome can be found here: https://www.genome.wisc.edu/sequencing/updating.htm

Where:
-f is the path to the genome fasta file
-m is the maximum number of mismatches that are allowed when determining a consensus sequence
-k is the length of the consensus sequence
-w is the length of the genome scanned near a GC-skew minimum
