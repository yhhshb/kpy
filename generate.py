#!/usr/bin/python3

import sys
import gzip
import random
import fastx

def uniform_main(args):
    random.seed(args.s)
    if args.o: fd = open(args.o, "w") 
    else: fd = sys.stdout
    name = 0
    for _ in range(0, args.n-1):
        seq = ''.join(random.choice('ACTG') for _ in range(args.l))
        fastx.fasta_write(fd, str(name), seq, False)
        name += 1
    if args.n >= 1:
        seq = ''.join(random.choice('ACTG') for _ in range(args.l))
        fastx.fasta_write(fd, str(name), seq, True)
    fd.flush()
    fd.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest = "command")

    parser_uniform = subparsers.add_parser("uniform", help = "Generate random uniform sequences")
    parser_uniform.add_argument("-l", help = "length of each sequence", type=int, required = True)
    parser_uniform.add_argument("-n", help = "number of sequences", type=int, required=True)
    parser_uniform.add_argument("-o", help = "fasta output [stdout]", type=str, required=False)
    parser_uniform.add_argument("-s", help = "random seed [42]", type=int, default=42)

    args = parser.parse_args(sys.argv)
    if (args.command == "uniform"): uniform_main(args)
    else: parser.print_help(sys.stderr)
    