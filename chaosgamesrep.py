import collections
import os, sys, time, re
import numpy as np
import argparse
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

    gKC_exec = "/home/amorris/BioInf/CGR/bin/genKmerCount"
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

    total_kmers = np.sum([float(c) for k, c in ksa])

    for k, c in ksa:
        freq_array[k] = float(c)/total_kmers

    return freq_array

def chaos_game_representation(probabilities, k, l):
    """ Where:
    A is upper left
    T is lower right
    G is upper right
    C is lower left
    """
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

def gen_plot(fasta_file, chaos_kn, k, l):
    a_size = len(chaos_kn)
    fig, ax = plt.subplots()
    im = ax.imshow(chaos_kn, interpolation='nearest', cmap=cm.gray_r)

    # ax.tick_params(axis='both',
    # which='both',
    # bottom=False,
    # top=False,
    # labelbottom=False,
    # right=False,
    # left=False,
    # labelleft=False)

    # ax.text(-(0.05*a_size), -(0.01*a_size), "A", fontsize=12)
    # ax.text(1.01*a_size, -(0.01*a_size), "G", fontsize=12)
    # ax.text(-(0.05*a_size), 1.05*a_size, "C", fontsize=12)
    # ax.text(1.01*a_size, 1.05*a_size, "T", fontsize=12)

    ax.set_title(f"{path.basename(path.splitext(fasta_file)[0])} k={k} log({l})")
    fig.tight_layout()
    plt.savefig(path.join(PATH, f"cgr_k{k}_log{np.round(l, 3)}.png"))
    # plt.show()

def is_file(filename):
    """ Checks if a path is a file """

    if not path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return path.abspath(path.realpath(path.expanduser(filename)))

def is_dir(direname):
    """ Checks if a path is a directory """

    if not path.isdir(direname):
        msg = "{0} is not a directory".format(direname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return path.abspath(path.realpath(path.expanduser(direname)))

def check_kmer_range(kmer_args):
    """ Checks given kmer args """

    try:
        kmer_args=[int(i) for i in kmer_args]
    except:
        raise(f"-k only takes integers! argerror: {i}")

    if len(kmer_args) == 1:
        print(f"Generating CGRs using the kmer size {kmer_args[0]}.")
        return kmer_args[0]
    elif len(kmer_args) == 2:
        if kmer_args[0] < kmer_args[1]:
            print(f"Kmer range must be i to j where i<j.")
            sys.exit(1)
        else:
            print(f"Generating CGRs using the kmer range {kmer_args[0]}-{kmer_args[1]}.")
            return kmer_args
    elif len(kmer_args) < 2:
        print(f"Too many arguments provided to -k. Kmer range must be i to j where i<j.")
        sys.exit(1)
    else:
        return [3]

def parse_args(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('script', type=path.abspath, action='store', help=argparse.SUPPRESS)
    parser.add_argument('fasta', type=is_file, action='store',
                        help='A FASTA format sequence file to transform into a chaos games represetation.')

    parser.add_argument('-k', nargs="*", type=int, default=[3, 13], action='store',
                        help='Kmer size for generating a chaos games representation. If two integers are provided, it carries out multiple games using kmer sizes over the given range. Default=3. \
                        Default=3.')
    parser.add_argument('-l', type=int, default=10, action='store',
                        help='A log value to transform data within the frequency matrix.')

    args = parser.parse_args(argv)
    return args

if __name__ == "__main__":

    args = parse_args(sys.argv)

    PATH = os.getcwd()

    start = time.time()

    logspace = np.logspace(1, 0.1, num=args.k[1]-args.k[0]+1, endpoint=True, base=args.l)

    for i, k in enumerate(range(args.k[0], args.k[1]+1)):
        log = logspace[i]/args.l

        print(f"k={k} log={log}")

        kc_handle = f"{args.fasta}.k{k}"

        if not path.isfile(kc_handle):
            kmer_count_file = kmerise(args.fasta, k, PATH)
        else:
            print(f"Kmer count file found at {kc_handle}")
            kmer_count_file = kc_handle

        freq_vector = read_kmer_array(kmer_count_file, k)
        print(f"Time elapsed: {time.time()-start}")

        chaos_kn = chaos_game_representation(freq_vector, k, log)
        gen_plot(args.fasta, chaos_kn, k, log)
