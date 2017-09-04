# Origin-of-Replication-Finder
A Bioinformatic Algorithm to Find the Origin of Replication on Bacterial Genomes

Because of the mechanics of bacterial division, GC-skew reaches a minimum around the origin of replication. This algorithm scans for consensus sequences around this region.

An example E coli genome can be found here: https://www.genome.wisc.edu/sequencing/updating.htm. This can be run from the terminal like so as './main.py -f E_coli_genome.fas -m 1 -k 9 -w 1000'.

Where:
* f is the path to the genome fasta file
* m is the maximum number of mismatches that are allowed when determining a consensus sequence
* k is the length of the consensus sequence
* w is the length of the genome scanned near a GC-skew minimum

