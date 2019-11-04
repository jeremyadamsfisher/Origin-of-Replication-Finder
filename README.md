# Origin of Replication Finder
### A Bioinformatic Algorithm to Find the Origin of Replication on Bacterial Genomes

Based on Pavel Pazners MOOC

Check it out [here](https://github.com/jeremyadamsfisher/Origin-of-Replication-Finder/blob/master/Origin_of_Replication_Finder.ipynb).

Because of the molecular mechanics of bacterial division, GC-skew reaches a minimum around the origin of replication. This algorithm scans for consensus sequences around this region.

An example E coli genome can be found here: https://www.genome.wisc.edu/sequencing/updating.htm.

There are several settings that can be tweaked in the Jupyter notebook to modify the search parameters.
 - *E_coli_genome.fas* is the assumed path to the genome fasta file
 - *max_mismatches_allowed* is the maximum number of mismatches that are allowed when determining a consensus sequence 
 - *k_mer_length* is the length of the consensus sequence
 - *window_length* is the length of the genome scanned near a GC-skew minimum
