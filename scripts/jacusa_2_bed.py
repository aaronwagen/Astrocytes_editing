#!/usr/bin/env python3

# GY210827

''' jacusa_to_bed '''


from argparse import ArgumentParser, FileType
from collections import Counter
from copy import copy
from csv import DictReader
from functools import cached_property
from itertools import permutations
from signal import signal, SIGPIPE, SIG_DFL
from sys import stdin, stderr


signal(SIGPIPE, SIG_DFL)  # Gracefully handle downstream PIPE closure


# Globals #####################################################################



# Classes #####################################################################


class Variant():

    def __init__(self, **kwargs):
        # import all fields from the csv.DictReader
        self.__dict__.update(kwargs)
        self.start = int(self.start)
        self.end = int(self.end)
        self.score = float(self.score)
        self.filter = self.filter != '*'
        # convert self.bases* fields to dicts of ints, make sums
        self.reps_1 = 0
        self.counts_1 = Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0})
        self.reps_2 = 0
        self.counts_2 = Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0})
        for k, v in self.__dict__.items():
            if k.startswith('bases'):
                self.__dict__[k] = dict(
                    zip(('A', 'C', 'G', 'T'),
                        (int(c) for c in v.split(','))))
                if k.startswith('bases1'):
                    self.reps_1 += 1
                    self.counts_1.update(self.__dict__[k])
                else:  # k.startswith('bases2')
                    self.reps_2 += 1
                    self.counts_2.update(self.__dict__[k])

    @cached_property
    def depth_1(self):
        return sum(self.counts_1.values())

    @cached_property
    def depth_2(self):
        if self.reps_2:
            return sum(self.counts_2.values())
        return None

    @cached_property
    def min_depth(self):
        if self.reps_2:
            return min(self.depth_1, self.depth_2)
        return self.depth_1

    @cached_property
    def edit(self):
        counts_1 = copy(self.counts_1)
        del(counts_1[self.ref])
        if self.reps_2:
            counts_2 = copy(self.counts_2)
            del(counts_2[self.ref])
            return sorted(
                (counts_1.most_common(1)[0], counts_2.most_common(1)[0]),
                key=lambda x: x[1], reverse=True)[0][0]
        return counts_1.most_common(1)[0][0]

    @cached_property
    def edit_type(self):
        return f'{self.ref}{self.edit}'

    @cached_property
    def n_edited(self):
        return max(
            sum(
                self.__dict__[f'bases1{i}'][self.edit] > 0
                for i in range(1, self.reps_1 + 1)),
            sum(
                self.__dict__[f'bases2{i}'][self.edit] > 0
                for i in range(1, self.reps_2 + 1)))

    @cached_property
    def edit_fraction_1(self):
        return self.counts_1[self.edit] / self.depth_1

    @cached_property
    def edit_fraction_2(self):
        if self.reps_2:
            return self.counts_2[self.edit] / self.depth_2
        return None

    @cached_property
    def max_edit_fraction(self):
        if self.reps_2:
            return max(self.edit_fraction_1, self.edit_fraction_2)
        return self.edit_fraction_1

    @cached_property
    def min_edit_fraction(self):
        if self.reps_2:
            return min(self.edit_fraction_1, self.edit_fraction_2)
        return self.edit_fraction_1

    @cached_property
    def fold(self):
        if self.reps_2:
            try:
                return self.max_edit_fraction / self.min_edit_fraction
            except ZeroDivisionError:
                return float('inf')
        return None

    def is_snp(self, ratio=0.5):
        return self.min_edit_fraction > ratio

    def __repr__(self):
        ''' Return a BED6 representation of the variant '''
        if self.reps_2:
            return \
                f'{self.contig}\t' \
                f'{self.start}\t' \
                f'{self.end}\t' \
                f'{self.edit_type}_{self.counts_1[self.edit]}/{self.depth_1}_{self.counts_2[self.edit]}/{self.depth_2}\t' \
                f'{self.score * -1.0 if self.edit_fraction_1 > self.edit_fraction_2 else self.score}\t' \
                f'{self.strand}'
        return \
            f'{self.contig}\t' \
            f'{self.start}\t' \
            f'{self.end}\t' \
            f'{self.edit_type}_{self.counts_1[self.edit]}/{self.depth_1}\t' \
            f'{self.score}\t' \
            f'{self.strand}'


###############################################################################


if __name__ == '__main__':

    parser = ArgumentParser(
        description='jacusa_to_bed',
        epilog='Coverts and filters JACUSA output format to BED6')
    parser.add_argument(
        'input', nargs='*', type=FileType('r'), default=[stdin],
        help='input file(s) (use "-" or leave blank for stdin)')
    parser.add_argument(
        '-e', '--edit', type=str,
        choices=[''.join(p) for p in permutations('ACGT', 2)],
        help='base transition to extract')
    parser.add_argument(
        '-z', '--z_score', type=float, default=1.960,
        help='threshold (absolute) Z-score (default %(default)s)')
    parser.add_argument(
        '-f', '--fold', type=float, default=2,
        help='threshold fold difference required between groups \
        (default %(default)s)')
    parser.add_argument(
        '-d', '--depth', type=int, default=10,
        help='minimum depth (summed for each group separately) required \
        (default %(default)s)')
    parser.add_argument(
        '-m', '--min_edit', type=float, default=0.01,
        help='minimum fraction of editing required in either group \
        (default %(default)s)')
    parser.add_argument(
        '-r', '--reps', type=int, default=1,
        help='minimum number of samples from either group that must display \
        editing (default %(default)s)')
    parser.add_argument(
        '-s', '--snp', type=int, default=0.5,
        help='treat any editing above this proportion as a snp and ignore \
        (default %(default)s)')
    parser.add_argument(
        '-q', '--quiet', action='store_true',
        help='silence reporting of the reasons for filtering')

    args = parser.parse_args()

    discards = {
        'JACUSA2 filter': 0, 'Z score threshold': 0, 'Potential SNP': 0,
        'Wrong editing type': 0, 'Insufficient depth': 0,
        'Insufficient editing': 0, 'Insufficient replicates with editing': 0,
        'Insufficient fold change': 0}
    passes = 0

    print(
        f'track '
        f'name=JACUSA2 '
        f'description="Stranded RNA modifications" '
        f'useScore=1')

    for f in args.input:
        _ = f.readline()  # skip program invocation line
        header = f.readline().lstrip('#').strip().split()
        reader = DictReader(f, header, dialect='excel-tab')
        for row in reader:
            v = Variant(**row)
            discard = False
            if v.filter:
                discards['JACUSA2 filter'] += 1
                continue
            if v.score <= args.z_score:
                discards['Z score threshold'] += 1
                continue
            if v.is_snp(args.snp):
                discards['Potential SNP'] += 1
                continue
            if args.edit and v.edit_type != args.edit:
                discards['Wrong editing type'] += 1
                continue
            if v.min_depth < args.depth:
                discards['Insufficient depth'] += 1
                continue
            if v.max_edit_fraction < args.min_edit:
                discards['Insufficient editing'] += 1
                continue
            if v.n_edited < args.reps:
                discards['Insufficient replicates with editing'] += 1
                continue
            if v.fold is not None and v.fold < args.fold:
                discards['Insufficient fold change'] += 1
                continue
            passes += 1
            print(v)

    if not args.quiet:
        print(
            f'{passes + sum(discards.values()):,d} sites inspected',
            file=stderr)
        for i, (k, v) in enumerate(discards.items(), 1):
            print(f'{".. " * i}{v:,d}: {k}', file=stderr)
        print(f'{passes:,d} sites passed filtering', file=stderr)

