#!/usr/bin/python3

"""
This script takes a FASTA file in input (-i option or stdin) and prints the randomly modified records.
"""

import sys
import gzip
import random
from enum import Enum
import fastx

class MType(Enum):
    SUBSTITUTION = 0
    INSERTION = 1
    DELETION = 3
    IDENTITY = 4

def mutate(seq: str, mut_rate: float, indel_fraction: float, ext_prob: float) -> str:
    '''Mutate the given sequence.
    '''
    nseq = list()
    indel_type = MType.IDENTITY
    indel_len = 0
    for base in seq:
        c = fastx.enc_nuc_table[ord(base)]
        #sys.stderr.write("c: {} -> ".format(c))
        if (indel_len == 0 and c < 4 and random.random() < mut_rate):
            if (random.random() >= indel_fraction):#substitution
                #sys.stderr.write("[muttype]: substitution -> ")
                r = random.random()
                c = (c + int(r * 3.0 + 1)) & 3
                indel_type = MType.SUBSTITUTION
            else:#indel
                #sys.stderr.write("[muttype]: ")
                indel_len = 1
                while(random.random() < ext_prob):#bernoulli trials
                    indel_len += 1
                if (random.random() < 0.5):
                    #sys.stderr.write("deletion -> ")
                    indel_type = MType.DELETION
                else:
                    #sys.stderr.write("insertion -> ")
                    indel_type = MType.INSERTION
        if (indel_type == MType.INSERTION):
            #sys.stderr.write("[apply] insertion l = {}".format(indel_len))
            while(indel_len > 0):
                nseq.append(fastx.dec_nuc_table[random.randrange(5)])
                indel_len -= 1
        elif (indel_type == MType.DELETION):
            #sys.stderr.write("[apply] deletion l = {}".format(indel_len))
            indel_len -= 1
        elif (indel_type == MType.IDENTITY or indel_type == MType.SUBSTITUTION):
            #sys.stderr.write("[apply] substitution or identity c = {}".format(c))
            nseq.append(fastx.dec_nuc_table[c])
        if(indel_len == 0):
            #sys.stderr.write(" -> [reset]")
            indel_type = MType.IDENTITY
        #sys.stderr.write("\n")
    return ''.join(nseq)

def fixed_main(args):
    random.seed(args.s)
    if args.i == None: fd = sys.stdin
    elif fastx.is_gzipped(args.i): fd = gzip.open(args.i, "rt")
    else: fd = open(args.i, "r")

    if (args.o == None): fo = sys.stdout
    else: fo = open(args.o, "w")
    prev = [None, None]
    for name, seq, _ in fastx.read(fd):
        if prev != [None, None]: 
            sys.stderr.write("{}|end\n".format(prev[0]))
            sys.stderr.write("{}\n".format(prev[1][-20:]))
            fastx.fasta_write(fo, prev[0], mutate(prev[1], args.m, args.d, args.x), False)#write buffer
        prev[0] = name
        prev[1] = seq
    if prev != [None, None]: 
        sys.stderr.write("{}|end\n".format(prev[0]))
        sys.stderr.write("{}\n".format(prev[1][-20:]))
        fastx.fasta_write(fo, prev[0], mutate(prev[1], args.m, args.d, args.x), True)
    fo.flush()
    fo.close()
    fd.close()

def random_main(args):
    random.seed(args.s)
    if args.i == None: fd = sys.stdin
    elif fastx.is_gzipped(args.i): fd = gzip.open(args.i, "rt")
    else: fd = open(args.i, "r")

    if (args.o == None): fo = sys.stdout
    else: fo = open(args.o, "w")
    if (args.log):
        fl = open(args.log, "w")
        fl.write("mutation rate,InDel fraction,ext. prob.\n")
    else:
        fl = None
    mrrange = (min(args.m), max(args.m))
    idrange = (min(args.d), max(args.d))
    eprange = (min(args.x), max(args.x))
    prev = [None, None]
    for id, seq, _ in fastx.read(fd):
        if (mrrange[0] == mrrange[1]): mutation_rate = mrrange[0]
        else: mutation_rate = random.uniform(mrrange[0], mrrange[1])
        if (idrange[0] == idrange[1]): indel_fraction = idrange[0]
        else: indel_fraction = random.uniform(idrange[0], idrange[1])
        if (eprange[0] == eprange[1]): extension_prob = eprange[0]
        else: extension_prob = random.uniform(eprange[0], eprange[1])
        if prev != [None, None]: 
            fastx.fasta_write(mutate(fo, prev[0], prev[1], mutation_rate, indel_fraction, extension_prob))
            if (fl): fl.write(str(mutation_rate)+'\t'+str(indel_fraction)+'\t'+str(eprange)+'\n')
        prev[0] = id
        prev[1] = seq
    if prev != [None, None]:
        if (mrrange[0] == mrrange[1]): mutation_rate = mrrange[0]
        else: mutation_rate = random.uniform(mrrange[0], mrrange[1])
        if (idrange[0] == idrange[1]): indel_fraction = idrange[0]
        else: indel_fraction = random.uniform(idrange[0], idrange[1])
        if (eprange[0] == eprange[1]): extension_prob = eprange[0]
        else: extension_prob = random.uniform(eprange[0], eprange[1])
        fastx.fasta_write(mutate(fo, prev[0], prev[1], mutation_rate, indel_fraction, extension_prob))
        if (fl): fl.write(str(mutation_rate)+'\t'+str(indel_fraction)+'\t'+str(eprange)+'\n')
    fd.close()
    fo.flush()
    fo.close()

def main(args):
    if (args.command == "fixed"): fixed_main(args)
    elif (args.command == "random"): random_main(args)
    else: parser.print_help(sys.stderr)
    
def setup_parser():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest = "command")

    parser_fixed = subparsers.add_parser("fixed", help="mutate all sequences in the FASTX file with the given parameters")
    parser_fixed.add_argument("-i", help="input FASTX file to modify [stdin as fasta]", type=str)
    parser_fixed.add_argument("-o", help="output FASTX file [stdout as fasta]", type=str)
    parser_fixed.add_argument("-m", help="mutation probability (substitutions + indels) [0.05]", type=float, default=0.05)
    parser_fixed.add_argument("-d", help="indel fraction of the mutations [0.15]", type=float, default=0.15)
    parser_fixed.add_argument("-x", help="probability that an indel is extended [0.3]", type=float, default=0.3)
    parser_fixed.add_argument("-s", help="random seed [42]", type=int, default=42)

    parser_random = subparsers.add_parser("random", help="mutate the sequences in the FASTX file with random mutation rates")
    parser_random.add_argument("-i", help="input FASTX file to modify [stdin as fasta]", nargs='?')
    parser_random.add_argument("-o", help="output FASTX file [stdout as fasta]", nargs='?')
    parser_random.add_argument("-m", help="range of mutation probabilities [[0.0, 0.7]]", nargs='+', type=float, default=[0.0, 0.7])
    parser_random.add_argument("-d", help="range of indel fractions of the mutations [[0.15]]", nargs='+', type=float, default=[0.15])
    parser_random.add_argument("-x", help="range of probabilities that an indel is extended [[0.3]]", nargs='+', type=float, default=[0.3])
    parser_random.add_argument("-s", help="random seed [42]", type=int, default=42)
    parser_random.add_argument("-l", "--log", help="log file in tsv format containing the parameters for each sequence", required=False)

    return parser

if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args(sys.argv)
    main(args)