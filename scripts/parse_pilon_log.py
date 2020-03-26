#!/usr/bin/env python3

import sys
from pathlib import Path
import re

reads_pattern = re.compile(
    r'Total Reads: (?P<reads>\d+), '
    r'Coverage: (?P<coverage>\d+), '
    r'minDepth: (?P<mindepth>\d+)')
confirmed_pattern = re.compile(
    r'Confirmed (?P<confirmed>\d+) of (?P<bases>\d+) bases')
correct_pattern = re.compile(
    r'Corrected (?P<snps>\d+) snps; '
    r'(?P<ambiguous>\d+) ambiguous bases; '
    r'corrected (?P<small_ins>\d+) small insertions '
    r'totaling (?P<small_ins_bases>\d+) bases, '
    r'(?P<small_del>\d+) small deletions '
    r'totaling (?P<small_del_bases>\d+) bases')

stats_dict = {
    'contigs': 0,
    'confirmed_bases': 0,
    'total_bases': 0,
    'reads': 0,
    'snps': 0,
    'ins': 0,
    'ins_bases': 0,
    'del': 0,
    'del_bases': 0
}

for fname in snakemake.input:
    with open(fname) as f:
        for line in f:
            reads_match = re.match(reads_pattern, line)
            if reads_match:
                reads = int(reads_match.group('reads'))
                coverage = int(reads_match.group('coverage'))
                mindepth = int(reads_match.group('mindepth'))
                stats_dict['reads'] += reads
                stats_dict['contigs'] += 1
                continue

            confirmed_match = re.match(confirmed_pattern, line)
            if confirmed_match:
                confirmed = int(confirmed_match.group('confirmed'))
                total_bases = int(confirmed_match.group('bases'))
                stats_dict['confirmed_bases'] += confirmed
                stats_dict['total_bases'] += total_bases
                continue

            correct_match = re.match(correct_pattern, line)
            if correct_match:
                snps = int(correct_match.group('snps'))
                ins = int(correct_match.group('small_ins'))
                ins_bases = int(correct_match.group('small_ins_bases'))
                dels = int(correct_match.group('small_del'))
                dels_bases = int(correct_match.group('small_del_bases'))
                stats_dict['snps'] += snps
                stats_dict['ins'] += ins
                stats_dict['ins_bases'] += ins_bases
                stats_dict['del'] += dels
                stats_dict['del_bases'] += dels_bases
                continue

with open(snakemake.output[0], 'w') as of:
    for key, value in stats_dict.items():
        print('{key}\t{value}'.format(key=key, value=value), file=of)
