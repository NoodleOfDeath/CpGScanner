#!/usr/bin/python3

import math
import multiprocessing
import multiprocessing.pool
import random
import string

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class NDPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


class StringRepresentable:

    def __str__(self):
        pass

    def __repr__(self):
        return str(self)


class Chunk(StringRepresentable):

    def __init__(self, seq, score, length):
        self.seq = seq
        self.score = score
        self.length = length

    def __str__(self):
        return "({0}, {1}, {2})".format(self.seq, self.reduce(), self.length)

    def reduce(self):
        return self.score / self.length


class Island(StringRepresentable):

    def __init__(self, seq, start, length):
        self.seq = seq
        self.start = start
        self.length = length

    def __str__(self):
        return "({0}, {1}, {2})".format(self.seq, self.start, self.length)


# Finds CpG islands based on passed criteria.
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
    dn = 0
    for i in range(0, len(chunks)):
        chunk = chunks[i]
        if chunk.reduce() < threshold or i == len(chunks) - 1:
            if dn >= min_length:
                length = offset - n
                island = Island(seq[n:n+length], n, offset - n)
                islands.append(island)
            n = -1
            dn = 0
        else:
            #print(chunk.seq)
            if n < 0:
               n = offset
            else:
               dn += 1
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
        with NDPool(threads) as pool:
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
            'threshold': .50,
            'min-length': 2,
        })
    print(islands)