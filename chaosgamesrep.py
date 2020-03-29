import collections
import os, sys, time, re
import numpy as np
from os import path, system
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
from Bio import SeqIO

import pylab
import math

def kmerise(fasta, k, PATH):
    """ Takes a sequence file in FASTA format and digests it into a counted kmer spectrum using genKmerCount.
    Outputs a tab seperated kmer count file.

    Source arguments used in this function:
    script = the path to this script, used to locate the genKmerCount executable
    k = kmer size used to digest the fasta file
    m = maximum kmer count to return kmers
    """

    gKC_exec = "/home/amorris/software/TreeMer/bin/genKmerCount"
    ks_outfile = path.join(PATH, f"{path.basename(fasta)}.k{k}")
    gKC_argline = f"{gKC_exec} {fasta} {k} {0} > {ks_outfile}"

    ## Spawn subprocess
    system(gKC_argline)

    return ks_outfile

def read_kmer_array(kmer_count_file, args):

    freq_array = OrderedDict()
    ksa = []

    with open(kmer_count_file, "r") as f:
        argc = f.readline()
        header = f.readline()

        for line in f.readlines():
            kmer, count, _ = re.split("[\t|\n]", line)
            if "N" in kmer:
                continue
            # kmer_array[kmer] = int(count)
            ksa.append([kmer, count])

    total_kmers = np.sum([int(c) for k, c in ksa])

    for k, c in ksa:
        freq_array[k] = float(c)/total_kmers

    return freq_array

def chaos_game_representation(probabilities, k, l):
    array_size = int(math.sqrt(4**k))
    chaos = []
    for i in range(array_size):
        chaos.append([0]*array_size)

    maxx = array_size
    maxy = array_size
    posx = 1
    posy = 1
    for key, value in probabilities.items():
        for char in key:
            if char == "T":
                posx += maxx / 2
            elif char == "C":
                posy += maxy / 2
            elif char == "G":
                 posx += maxx / 2
                 posy += maxy / 2
            maxx = maxx / 2
            maxy /= 2

        chaos[int(posy)-1][int(posx)-1] = value**l
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
    return chaos

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("USAGE: chaosgamesrep.py <fasta> <k> <logscale>")
        sys.exit(1)

    PATH = os.getcwd()

    fasta = sys.argv[1]
    k = int(sys.argv[2])
    l = float(sys.argv[3])

    start = time.time()
    kmer_count_file = kmerise(fasta, k, PATH)
    # kmer_count_file = "/home/amorris/CGR/byzA.cds.k13"
    freq_vector = read_kmer_array(kmer_count_file, k)
    print(f"Time elapsed: {time.time()-start}")

    chaos_kn = chaos_game_representation(freq_vector, k, l)
    pylab.title(f"CGR k={k} log({l})")
    pylab.imshow(chaos_kn, interpolation='nearest', cmap=cm.gray_r)
    pylab.savefig(os.path.join(PATH, f"chaos{k}_log{l}.png"))
    pylab.show()
