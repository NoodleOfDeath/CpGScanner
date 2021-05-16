#!/usr/bin/python3

import math
import multiprocessing
import multiprocessing.pool
import random
import string

# Define a no-daemon process (NoDaemonProcess) class so we can have processes spawned by a pool and will spawn their own subprocesses
class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


# Define a no-daemon pool (NoDaemonPool) class that spawns no-daemon processes which will further spawn their own subprocesses
class NoDaemonPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


# Define a string representable (StringRepresentable) interface for debug purposes
class StringRepresentable:

    def __str__(self):
        pass

    def __repr__(self):
        return str(self)


# Define Chunk class for parallel processing of the genome
class Chunk(StringRepresentable):

    def __init__(self, seq, score, length):
        self.seq = seq
        self.score = score
        self.length = length

    def __str__(self):
        return "({0}, {1}, {2})".format(self.seq, self.reduce(), self.length)

    def reduce(self):
        return self.score / self.length


# Define a CpG island (Island) class that represents a CpG site in a genomic sequence, in the form of a substring with an index and length
class Island(StringRepresentable):

    def __init__(self, seq, index, length):
        self.seq = seq
        self.index = index
        self.length = length

    def __str__(self):
        return "({0}, {1}, {2})".format(self.seq, self.index, self.length)


# Finds CpG islands in a passed genomic sequence based on passed criteria.
def find_islands(seq, opt = {}):
    if not 'threshold' in opt:
        print('Missing threshold option')
        return
    if not 'min-length' in opt:
        print('Missing min-length option')
        return
    threshold = opt['threshold']
    min_length = opt['min-length']
    chunks = seek(seq, opt)
    print("Sliced %ld character sequence into %ld chunks" % (len(seq), len(chunks)))
    #print(chunks)
    # Search for islands of min_length that meet the threshold
    islands = []
    offset = 0
    n = -1
    for i in range(0, len(chunks)):
        chunk = chunks[i]
        if chunk.reduce() < threshold or i == len(chunks) - 1:
            if n >= 0 and offset - n >= min_length:
                length = offset - n
                island = Island(seq[n:n+length], n, offset - n)
                islands.append(island)
            n = -1
        else:
            #print(chunk.seq)
            if n < 0:
               n = offset
        offset += chunk.length
    print("Found %ld CpG islands matching the criteria" % len(islands))
    return islands

# Subroutine for breaking genetic sequence into chunks to be processed in parallel.
def seek(seq, opt = {}):
    threads = 8
    if 'threads' in opt:
        threads = opt['threads']
    if not 'chunk' in opt:
        print('Missing chunk option')
        return
    chunk_size = opt['chunk']
    #print(seq, len(seq))
    if len(seq) > chunk_size:
        mid = math.floor(len(seq)/2)
        with NoDaemonPool(threads) as pool:
            results = []
            # Process chunks in parallel
            products = pool.starmap(seek, [(seq[0:mid], opt), (seq[mid:len(seq)], opt)])
            if not type(products[0]) is Chunk:
                results += products[0]
                results += products[1]
            else:
                results += products
            return results
    else:
        score = seq.count('C') + seq.count('G')
        return Chunk(seq, score, len(seq))

def gen_seq(length):
    return ''.join(random.choice(['C', 'G', 'A', 'T']) for i in range(length))

if __name__ == '__main__':
    islands = find_islands(
        gen_seq(256),
        { 
            'threads': 8,
            'chunk': 4,
            'threshold': 0.6,
            'min-length': 8,
        })
    print(islands)