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

    def __next__(self):
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