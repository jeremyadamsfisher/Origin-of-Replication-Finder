import regex
from pandas import DataFrame
from collections import defaultdict
        
## CRITICAL CODE -- CRUCIAL METHODS FOR THE IMPLEMENTATION OF THE ORIC FINDER ALGORITHM ##

def calculate_minimum_skew_locations(_genome):
    '''
    Calculates the GC skew curve of the genome and returns
    its absolute minimums
    '''
    
    nucleotide_to_delta_skew = {'A': 0, 'C': -1, 'T': 0, 'G': 1, 'N': 0}
    
    df = DataFrame({'Nucleotide': [n for n in _genome],
                    'Delta_Skew': [nucleotide_to_delta_skew[n] for n in _genome]})
    df['Skew'] = df['Delta_Skew'].cumsum()
    minimum_skew = df['Skew'].min()
    minimum_skew_locations = [int(l) for l in df[df['Skew'] == minimum_skew].index if df.loc[l]['Nucleotide'] == 'C']
    
    return minimum_skew_locations

def window_centered_around(_center, _window_length, _genome):
    '''
    Returns the polymer centered around a specified point, of a
    specified length, within the given genome
    '''
    
    return _genome[int(_center - _window_length / 2):int(_center + (_window_length / 2))]

def most_frequent_kmers(_sequence, _k_mer_length, _max_mismatches):
    '''
    Returns the kmers with the greatest number of hits and the
    number of those hits
    '''
    
    def sequence_neighborhood(pattern, max_mismatches_allowed):
        '''
        Returns a 'neighborhood' of sequences where there are a defined
        number of mismatches with the seed sequence
        '''
    
        substitutes = {'A': ['C', 'T', 'G'],
                       'T': ['A', 'C', 'G'],
                       'G': ['A', 'T', 'C'],
                       'C': ['A', 'T', 'G'],
                       'N': ['C', 'A', 'T', 'G']}
    
        def sequences_with_one_mismatch(pattern):

            sequences = set([pattern])

            for i, nucleotide in enumerate(pattern):
                for substitute_nucleotide in substitutes[nucleotide]:
                    new_neighbor = pattern[:i] + substitute_nucleotide + pattern[i + 1:]
                    sequences.add(new_neighbor)

            return list(sequences)

        neighborhood = set([pattern])

        for i in range(max_mismatches_allowed):
            new_neighbors = [sequences_with_one_mismatch(neighbor) for neighbor in neighborhood]
            neighborhood.update(flatten(new_neighbors))

        return neighborhood
        
    def num_approx_matches(_in_sequence, _of_kmer, max_mismatches_allowed):
        '''
        Compares two sequence and returns the number of times there is an
        aproximate match, defined by having at most a defined number of
        mismatches
        '''
        search_expression = '(%s){s<=%s}' % (_of_kmer, max_mismatches_allowed)
        occurrences = regex.findall(search_expression, _in_sequence)
        return len(occurrences)
        
        
        
    
    possible_kmers = set([])

    for sequence in Subsequences(_sequence, k):
        possible_kmers.update(sequence_neighborhood(sequence, max_mismatches_allowed))

    k_mer_hits = defaultdict(list)
    
    for kmer in possible_kmers:
        forward_hits = num_approx_matches(_sequence, kmer, _max_mismatches)
        reverse_complement_hits = num_approx_matches(_sequence, reverse_complement(kmer), _max_mismatches)
        k_mer_hits[forward_hits + reverse_complement_hits].append(kmer)
        
    number_of_top_hits = max(k_mer_hits.keys())
    top_hits = k_mer_hits[number_of_top_hits]
    
    return number_of_top_hits, top_hits
    
## SUPPORT CODE -- USEFUL METHODS THAT ARE NOT CRITICAL FOR THE IMPLEMENTATION ##

class Subsequences:
    '''
    Subsequences is an iterator that iterates over all
    possible polymers of a defined length within a sequence
    '''
    def __init__(self, _sequence, _kmer_length):
        self.sequence = _sequence
        self.k = _kmer_length
        self.pos = 0

    def __iter__(self):
        return self

    def __len__(self):
        return len(self.sequence) - self.k + 1

    def __getitem__(self, item):
        return self.sequence[item:item+self.k]

    def next(self):
        if self.pos < len(self.sequence) - self.k:
            self.pos += 1
            return self.sequence[self.pos:self.pos+self.k]
        else:
            raise StopIteration
            
def flatten(l):
    '''
    Flatten a nestled list
    '''
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out


def reverse_complement(sequence):
    '''
    Returns the reverse complement of the sequence
    '''
    complement = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T', 'N': 'N'}
    try:
        return ''.join(complement[bp] for bp in reversed(sequence))
    except KeyError:
        raise Exception('Attempted to find the complement of a base pair that is not A, C, G, T or N.')
