from functions import *

def main(max_mismatches_allowed, window_length, k_mer_length, genome_file):
    
    with open(genome_file, 'r') as f:

        # Read genome, typically as a FASTA file
        genome = ''.join(line.strip().upper() for line in f if not line.startswith('>'))
        
        # Calculate the absolute minimum skew locations
        for min_skew_loc in minimum_skew_locations(genome):
            
            # Search for the polymers with the greatest number of (aproximate) matches
            # within a window centered around this minimum skew location
            window = window_centered_around(min_skew_loc, window_length, genome)
            number_of_top_hits, top_hits = most_frequent_kmers(window, k_mer_length, max_mismatches_allowed):

            # Print the number of hits and the relevant polymer to the console
            if(number_of_top_hits > 0):
                print('Found the following polymers centered around {}:'.format(min_skew_loc))
                for hit in top_hits:
                    print('\t* {} [{} hits]'.format(hit, number_of_top_hits))