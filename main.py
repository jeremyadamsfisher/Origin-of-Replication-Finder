#! /usr/bin/python

import regex
from pandas import DataFrame
from collections import defaultdict
from optparse import OptionParser

nucleotide_to_delta_skew = {'A': 0, 'C': -1, 'T': 0, 'G': 1, 'N': 0}

complement = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T', 'N': 'N'}

substitutes = {'A': ['C', 'T', 'G'],
               'T': ['A', 'C', 'G'],
               'G': ['A', 'T', 'C'],
               'C': ['A', 'T', 'G'],
               'N': ['C', 'A', 'T', 'G']}


class Subsequences:
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
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out


def reverse_complement(sequence):
    try:
        return ''.join(complement[bp] for bp in reversed(sequence))
    except KeyError:
        raise Exception('Attempted to find the complement of a base pair that is not A, C, G, T or N.')


def sequence_neighborhood(pattern, max_mismatches_allowed):
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


def num_approx_matches(_in, _of, max_mismatches_allowed):
    search_expression = '(%s){s<=%s}' % (_of, max_mismatches_allowed)
    occurrences = regex.findall(search_expression, _in)
    return len(occurrences)


def main(max_mismatches_allowed, window_length, k, genome_file):
    with open(genome_file, 'r') as f:
        genome = f.read()

        df = DataFrame({'Nucleotide': [n for n in genome],
                        'Delta_Skew': [nucleotide_to_delta_skew[n] for n in genome]})
        df['Skew'] = df['Delta_Skew'].cumsum()

        minimum_skew = df['Skew'].min()
        minimum_skew_locations = [int(l) for l in df[df['Skew'] == minimum_skew].index if df.loc[l]['Nucleotide'] == 'C']

        for m in minimum_skew_locations:
            print('Calculating frequent 9-mers around location {}:'.format(m + 1))
            window = genome[m - (window_length / 2):m + (window_length / 2)]

            possible_kmers = set([])

            for sequence in Subsequences(window, k):
                possible_kmers.update(sequence_neighborhood(sequence, max_mismatches_allowed))

            k_mer_hits = defaultdict(list)
            for kmer in possible_kmers:
                forward_hits = num_approx_matches(window, kmer, max_mismatches_allowed)
                reverse_complement_hits = num_approx_matches(window, reverse_complement(kmer), max_mismatches_allowed)
                k_mer_hits[forward_hits + reverse_complement_hits].append(kmer)

            number_of_top_hits = max(k_mer_hits.keys())
            top_hits = k_mer_hits[number_of_top_hits]

            for hit in top_hits:
                print('\t* {} [{} hits]'.format(hit, number_of_top_hits))

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-f', '--file', dest='genome_file',
                      help='Path to the genome file', metavar='FILE')
    parser.add_option('-k', '--kmer-length', dest='kmer_length',
                      help='Length of the consensus desired', metavar='KMER')
    parser.add_option('-m', '--mismatches-allowed', dest='num_mismatches',
                      help='The number of mutations allowed in the consensus sequence', metavar='MIS')
    parser.add_option('-w', '--window-length', dest='win_len',
                      help='Length of the window around the origin of replication to search for consensus sequences', metavar='WIN')

    (options, args) = parser.parse_args()

    num_mismatches = int(options.num_mismatches)
    win_len = int(options.win_len)
    kmer_length = int(options.kmer_length)
    genome_file = options.genome_file

    main(num_mismatches, win_len, kmer_length, genome_file)
