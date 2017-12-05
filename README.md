# Origin of Replication Finder
### A Bioinformatic Algorithm to Find the Origin of Replication on Bacterial Genomes

Because of the molecular mechanics of bacterial division, GC-skew reaches a minimum around the origin of replication. This algorithm scans for consensus sequences around this region.

An example E coli genome can be found here: https://www.genome.wisc.edu/sequencing/updating.htm.

There are several settings that can be tweaked in the Jupyter notebook to modify the search parameters.
 - *f* is the path to the genome fasta file
 - *m* is the maximum number of mismatches that are allowed when determining a consensus sequence 
 - *k* is the length of the consensus sequence
 - *w* is the length of the genome scanned near a GC-skew minimum

This algorithm requires the *regex* module, which can be obtained by executing the following in the terminal:
~~~~shell
pip install regex
~~~~

