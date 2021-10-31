import builtins
import sys
import gzip
import math
import fastx

def smart_numeric_cast(s):
    def is_number(s: str):
        try:
            float(s)
            return True
        except ValueError:
            return False
    if is_number(s):
        n = float(s)
        if n.is_integer(): return int(s)#with large numbers casting float gives an approximation error, better to use the original string
        else: return n
    else:
        return s

def nibble2hex(msn: int, lsn: int):
    if msn >= 4 or lsn >= 4: raise ValueError("Bases must be 2bit encoded")
    n = msn * 4 + lsn
    return "{0:01X}".format(n)

def hex(km: str, enc_table: list) -> str:
    """Encode kmer as a 2bit-packed hex string"""
    en = [enc_table[ord(c)] for c in km]
    if len(en) % 2 == 1: en.append(0)
    l = len(en)
    res = list()
    print(en)
    for i in range(0, int(l / 2)): res.append(nibble2hex(en[i*2], en[i*2+1]))
    return "".join(res)

def set(seq: str, k: int, canonical: bool, table: set[str]) -> set[str]:
    l = len(seq)
    if l < k: return
    if table == None: table = builtins.set()
    for i in range(l - k + 1):
        kmer = seq[i:(i+k)]
        ok = True
        for c in kmer: 
            if not (c in fastx.base_for): ok = False
        if ok:
            if canonical:
                kmer_rev = kmer.translate(fastx.comp_tab)[::-1]
                if kmer > kmer_rev: kmer = kmer_rev
            table.add(kmer)
    return table

def count(seq: str, k: int, canonical: bool, table: dict) -> dict:
    l = len(seq)
    if l < k: return
    if table == None: table = dict()
    for i in range(l - k + 1):
        kmer = seq[i:(i+k)]
        ok = True
        for c in kmer: 
            if not (c in fastx.base_for): ok = False
        if ok:
            if canonical:
                kmer_rev = kmer.translate(fastx.comp_tab)[::-1]
                if kmer > kmer_rev: kmer = kmer_rev
            if kmer in table: table[kmer] += 1
            else: table[kmer] = 1
    return table

def diff(table1, table2, symmetric: bool = False) -> list:
    s1 = isinstance(table1, builtins.set)
    s2 = isinstance(table2, builtins.set)
    d1 = isinstance(table1, dict)
    d2 = isinstance(table2, dict)
    assert s1 or d1
    assert s2 or d2
    res = list()
    if s1: var1 = dict.fromkeys(table1, 1)
    else: var1 = table1
    if s2: var2 = dict.fromkeys(table2, 1)
    else: var2 = table2
    for key, val in var1.items():
        if key not in var2: 
            res.append(('i', key, val))
        else: 
            delta = val - var2[key]
            if delta: res.append(('i', key, delta))
    if symmetric:
        for key, val in var2.items():
            if key not in var1: 
                res.append(('j', key, val))
            else: 
                delta = val - var1[key]
                if delta: res.append(('j', key, delta))
    return res

def length(obj, sep=' ') -> int:
    import io
    if isinstance(obj, str):
        return len(obj)
    elif isinstance(obj, io.TextIOWrapper):
        kmer, _ = obj.readline().split(sep)
        obj.seek(0)
        return len(kmer)
    else:
        raise TypeError("Cannot handle object of type {}".format(type(obj)))

def minimizer(m: int, seq: str, hasher) -> str:
    l = len(seq)
    if l < m: return None
    minh = math.inf
    mini = None
    for i in range(l - m + 1):
        hval = hasher(seq[i:(i+m)].encode())
        if minh > hval: 
            minh = hval
            mini = i
    return seq[mini:(mini+m)], seq[:mini]+'|'+seq[mini+m:]

class Spectrum:
    
    def __init__(self):
        self.histogram = dict()
        
    def add(self, count: int):
        if count in self.histogram: self.histogram[count] += 1
        else: self.histogram[count] = 1

    def addFromDict(self, table: dict):
        for _, count in table.items():
            if count in self.histogram: self.histogram[count] += 1
            else: self.histogram[count] = 1

    def addFromFile(self, file, sep: str):
        if isinstance(file, str):
            file_handle = open(file, "r")
            close = True
        else: 
            file_handle = file
            close = False
        for line in file_handle:
            _, count = line.split(sep)
            count = smart_numeric_cast(count)
            if count in self.histogram: self.histogram[count] += 1
            else: self.histogram[count] = 1
        if close: file_handle.close()

    def empty(self):
        return not bool(self.histogram)

    def L0Norm(self):
        return sum(list(self.histogram.values()))

    def L1Norm(self):
        L1 = 0
        for k, v in self.histogram.items(): L1 += smart_numeric_cast(k) * smart_numeric_cast(v)
        return L1

    def entropy(self):
        #L1 = self.L1Norm()
        L0 = self.L0Norm()
        H0 = 0
        for _, column in self.histogram.items():
            p = column/L0
            H0 += -p*math.log(p,2)
        return H0

    def getMaxColumn(self):
        return max(list(self.histogram.values()))

    def getMaxCount(self):
        mcounts = list()
        colmax = self.getMaxColumn()
        for count, column in self.histogram.items():
            if column == colmax: mcounts.append(count)
        mcounts.sort()
        return mcounts[0]

    def getOptimalEpsilon(self):
        """Dimension a Bloom Filter for the cold items of the spectrum.

        :return: The best epsilon for a Bloom Filter storing cold items only in order to minimize the false positives coming from the most common item.
        """
        L0 = self.L0Norm()
        N = self.getMaxColumn()
        return (L0 - N) / N
        #if N == None:
        #    self.getMaxCount()
        #    N = self._L1light
        #L1 = self.L1Norm()
        #return N / (L1 - N)

    def getOptimizedEpsilon(self, c_csf: float):
        L0 = self.L0Norm()
        N = self.getMaxColumn()
        c_bf = 1/math.log(2)
        epsilon = c_bf/c_csf * ((L0 - N)/N) * math.log(math.e, 2)
        if epsilon == 0: epsilon = math.inf #raise RuntimeError("epsilon = 0 | L0 = {} | N = {} | c_csf = {}".format(L0, N, c_csf))
        return epsilon

    def removeCount(self, count: int):
        if count in self.histogram: 
            del self.histogram[count]

def set_main(args):
    table = builtins.set()
    for f in args.i:
        if f.endswith(".gz"): fi = gzip.open(f, "rt")
        else: fi = open(f, "r")
        for _, seq, _ in fastx.read(fi):
            set(seq, args.k, args.c, table)
        fi.close()
    if not args.i:
        for _, seq, _ in fastx.read(sys.stdin):
            set(seq, args.k, args.c, table)
    if (args.o): fo = open(args.o, "w")
    else: fo = sys.stdout
    for k in table: fo.write("{}\n".format(k))
    fo.close()

def count_main(args):
    table = dict()
    for f in args.i:
        if f.endswith(".gz"): fi = gzip.open(f, "rt")
        else: fi = open(f, "r")
        for _, seq, _ in fastx.read(fi):
            count(seq, args.k, args.c, table)
        fi.close()
    if not args.i:
        for _, seq, _ in fastx.read(sys.stdin):
            count(seq, args.k, args.c, table)
    if (args.o): fo = open(args.o, "w")
    else: fo = sys.stdout
    for k, v in table.items(): fo.write("{} {}\n".format(k, v))
    fo.close()

def diff_main(args):
    table1 = dict()
    table2 = dict()
    with open(args.i, "r") as th:
        for _, seq, _ in fastx.read(th):
            count(seq, args.k, args.c, table1)
    with open(args.j, "r") as th:
        for _, seq, _ in fastx.read(th):
            count(seq, args.k, args.c, table2)
    difference = diff(table1, table2, args.s)
    for s, km, delta in difference:
        sys.stdout.write("{},{},{}\n".format(s, km, delta))

def main(args):
    if (args.command == "set"): set_main(args)
    elif (args.command == "count"): count_main(args)
    elif (args.command == "diff"): diff_main(args)
    else: parser.print_help(sys.stderr)
    
def setup_parser():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("__default")
    subparsers = parser.add_subparsers(dest = "command")

    parser_set = subparsers.add_parser("set", help="Compute set of k-mers")
    parser_set.add_argument("-i", help="input FASTX file [stdin]", type=str, nargs='+', default=[])
    parser_set.add_argument("-o", help="output (one k-mer per line) [stdout]", type=str)
    parser_set.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_set.add_argument("-c", help="canonical k-mers", action="store_true")

    parser_count = subparsers.add_parser("count", help="Count k-mers")
    parser_count.add_argument("-i", help="input files (fasta or fastq) [stdin]", type=str, nargs='+', default=[])
    parser_count.add_argument("-o", help="output count table [stdout]", type=str)
    parser_count.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_count.add_argument("-c", help="canonical k-mers", action="store_true")

    parser_diff = subparsers.add_parser("diff", help="Compute k-mer difference between two fastx files")
    parser_diff.add_argument("-i", help="first fasta file", type=str, required=True)
    parser_diff.add_argument("-j", help="second fasta file", type=str, required=True)
    parser_diff.add_argument("-k", help="k-mer length", type=int, required=True)
    parser_diff.add_argument("-c", help="canonical k-mers", action="store_true")
    parser_diff.add_argument("-s", help="symmetric difference", action="store_true")

    return parser

if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args(sys.argv)
    main(args)

