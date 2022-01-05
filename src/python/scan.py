#!/usr/bin/env python3

import math
import multiprocessing
import multiprocessing.pool
import random
import re
import sys
import timeit

# Daemon processes are important for parallel processing and improving computational
# performance of this script.

class NoDaemonProcess(multiprocessing.Process):
    """ 
    Define a no-daemon process (NoDaemonProcess) class so we can have processes spawned 
    by a pool that can also spawn their own subprocesses
    """
    
    @property
    def daemon(self):
        return False
    
    @daemon.setter
    def daemon(self, value):
        pass

class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess

class NoDaemonPool(multiprocessing.pool.Pool):
    """
    Define a no-daemon pool (NoDaemonPool) class that spawns no-daemon processes 
    which will also spawn their own subprocesses
    """
    
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(NoDaemonPool, self).__init__(*args, **kwargs)


class StringRepresentable:
    """
    Define a string representable (StringRepresentable) interface for debug purposes
    """

    def __str__(self):
        pass

    def __repr__(self):
        return str(self)


class Chunk(StringRepresentable):
    """
    Define Chunk class for parallel processing of the genome
    """

    def __init__(self, seq, score, length):
        self.seq = seq
        self.score = score
        self.length = length

    def __str__(self):
        return "({0}, {1}, {2})".format(self.seq, self.reduce(), self.length)

    def reduce(self):
        """
        Returns the ratio of CG presence in this chunk
        """
        
        return self.score / self.length


class Island(StringRepresentable):
    """
    Define a CpG island (Island) class that represents a CpG site in a genomic sequence, 
    in the form of a substring with an index and length
    """

    def __init__(self, seq, index, length):
        self.seq = seq
        self.index = index
        self.length = length

    def __str__(self):
        return "({0}, {1}, {2})".format(self.seq, self.index, self.length)


def find_islands(seq, opt = {}):
    """
    Finds CpG islands in a passed genomic sequence based on passed criteria.
    """

    # Parse options and/or use defaults
    threads = 2
    if 'threads' in opt:
        threads = opt['threads']
    if not 'chunk' in opt:
        print('Missing chunk option')
        return
    if not 'threshold' in opt:
        print('Missing threshold option')
        return
    if not 'min-length' in opt:
        print('Missing min-length option')
        return

    chunk_size = opt['chunk']
    threshold = opt['threshold']
    min_length = opt['min-length']

    print("Searching %ld character sequence for CpG islands using %ld threads, threshold of %.2f, min length of %ld, and slicing into chunks of size %ld or less" % (len(seq), threads, threshold, min_length, chunk_size))
    
    # Breaks sequence into chunks of equal size based on options
    chunks = seek(seq, opt)

    print("Sliced %ld character sequence into %ld chunks" % (len(seq), len(chunks)))

    # Search for islands of `min_length` that meet the threshold from `chunks`
    islands = []
    offset = 0
    n = -1
    for i in range(len(chunks)):
        chunk = chunks[i]
        # If chunk does not meet threshold or is the last chunk, terminate the 
        # current island (if there is one) and append it to the return array.
        # Then reset `n` which is the offset of the next possible CpG island.
        if chunk.reduce() < threshold or i == len(chunks) - 1:
            # Only consider the island valid if it is greater in length than
            # `min_length`.
            if n >= 0 and offset - n >= min_length:
                length = offset - n
                island = Island(seq[n:n+length], n, offset - n)
                islands.append(island)
            n = -1
        else:
            # If `n < 0` then the index of a CpG island has not yet been set
            # indicating the start of a possible new CpG island.
            if n < 0:
               n = offset
        offset += chunk.length

    print("Found %ld CpG islands matching the criteria" % len(islands))
    return islands


def seek(seq, opt = {}):
    """
    Subroutine for breaking genomic sequence into chunks to be processed in parallel.
    """

    threads = 2
    if 'threads' in opt:
        threads = opt['threads']
    if not 'chunk' in opt:
        print('Missing chunk option')
        return

    chunk_size = opt['chunk']

    if len(seq) > chunk_size:
        mid = math.floor(len(seq)/2)
        with NoDaemonPool(threads) as pool:

            chunks = []
            # Process the chunks in parallel:
            # Recursively calls `seek` on the head and tail halves of `seq`
            r = pool.starmap(seek, [(seq[0:mid], opt), (seq[mid:len(seq)], opt)])

            # This is equivalent to serially calling `seek` like below:
            # r = (seek(seq[0:mid], opt), seek(seq[mid:len(seq)], opt))

            # If the return object is two arrays of chunks, append their contents,
            # separately so we don't end up with multidimensional arrays
            if not type(r[0]) is Chunk:
                chunks += r[0]
                chunks += r[1]
            else:
                chunks += r

            return chunks

    # TODO: is regex or count more performant than `count` method?
    score = len(re.findall('[CG]', seq)) # seq.count('C') + seq.count('G')
    return Chunk(seq, score, len(seq))


def gen_seq(length):
    """
    Generates a random genomic sequence of a specified length.
    """
    
    return ''.join(random.choice(['C', 'G', 'A', 'T']) for i in range(length))

if __name__ == '__main__':

    args = []
    threads = 2
    chunk = 4
    threshold = 0.6
    min_length = 8
    xnsize = 1024
    seq = ''

    i = 0
    argv = sys.argv[1:]

    # Parse out command line arguments into options
    while len(argv) > 0:
        arg = argv.pop(0)
        if arg == "-h" or arg == "--help":
            print(
                "usage: ./scan.py [-t <threads>] [-th <threshold>] [-c <chunk-size>] [-m <min-length>] [-n <random-sequence-length>] [genomic-sequence]"
            )
            exit()
        elif arg == "-t" or arg == "--threads":
            threads = int(argv.pop(0))
        elif arg == '-c' or arg == '--chunk':
            chunk = int(argv.pop(0))
        elif arg == '-th' or arg == '--threshold':
            threshold = float(argv.pop(0))
        elif arg == '-m' or arg == '--min-length':
            min_length = int(argv.pop(0))
        elif arg == '-n':
            xnsize = int(argv.pop(0))
        else:
            args.append(arg)

    if len(args) > 0:
        seq = args[0]

    if len(seq) == 0:
        start = timeit.default_timer()
        seq = gen_seq(xnsize)
        stop = timeit.default_timer()
        print("Took %.2f second(s) to generate %ld character sequence" % (stop - start, len(seq)))
    # Record start time
    start = timeit.default_timer()
    # Find islands in seq based on options passed by command line
    islands = find_islands(
        seq,
        { 
            'threads': threads,
            'chunk': chunk,
            'threshold': threshold,
            'min-length': min_length,
        })
    # Record stop time for debug purposes
    stop = timeit.default_timer()
    print("Took %.2f second(s) to find %ld CpG islands in %ld character sequence" % (stop - start, len(islands), len(seq)))
    print(islands)
