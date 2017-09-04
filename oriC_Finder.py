#! /usr/bin/python

from oriC_Finder import main
from optparse import OptionParser

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
