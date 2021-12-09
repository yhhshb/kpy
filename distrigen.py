import sys
import math
import random
import bisect
import functools

class ZipfGenerator: 

    def __init__(self, number_of_values: int, skew: float): 
        # Calculate Zeta values from 1 to n:
        tmp = [1. / (math.pow(float(i), skew)) for i in range(1, number_of_values+1)] 
        zeta = functools.reduce(lambda sums, x: sums + [sums[-1] + x], tmp, [0]) 

        # Store the translation map: 
        self.distMap = [x / zeta[-1] for x in zeta] 

    def next(self): 
        # Take a uniform 0-1 pseudo-random value:
        u = random.random()  

        # Translate the Zipf variable:
        return bisect.bisect(self.distMap, u) - 1 

def zipfian_main(args):
    assert args.ndraws >= 0
    assert args.nvals >= 0
    assert args.skew >= 0
    zipf = ZipfGenerator(args.nvals, args.skew)
    ln = 0
    while(ln < args.ndraws):
        sys.stdout.write("{} {}\n".format(ln, zipf.next()))
        ln += 1

def uniform_main(args):
    assert args.ndraws >= 0
    assert args.nvals >= 0
    for i in range(args.ndraws):
        sys.stdout.write("{} {}\n".format(i, random.randrange(0, args.nvals)))

def main(args):
    if (args.command == "zipfian"): return zipfian_main(args)
    elif (args.command == "uniform"): return uniform_main(args)
    # elif (args.command == "simulate"): return simulate_csf_main(args)
    else: sys.stderr.write("-h to list available subcommands\n")

def parser_init():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest="command")

    parser_zipfian = subparsers.add_parser("zipfian", help="Simulate data following a zipfian distribution")
    parser_zipfian.add_argument("-n", "--ndraws", help="numer of draws", type=int, required=True)
    parser_zipfian.add_argument("-v", "--nvals", help="number of values to draw from", type=int, required=True)
    parser_zipfian.add_argument("-s", "--skew", help="skew value of the array of values", type=float, required=True)

    parser_uniform = subparsers.add_parser("uniform", help="Simulate data following a uniform distribution")
    parser_uniform.add_argument("-n", "--ndraws", help="number of draws", type=int, required=True)
    parser_uniform.add_argument("-v", "--nvals", help="number of values", type=int, required=True)

    return parser

if __name__ == "__main__":
    parser = parser_init()
    args = parser.parse_args(sys.argv)
    main(args)