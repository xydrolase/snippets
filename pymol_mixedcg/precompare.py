#!/usr/bin/env python

import sys
import os
import re
from itertools import combinations

def find_meta(raw):
    return re.compile(r'^.+ - (\d+)$').search(raw[1]).group(1), raw[5]

if len(sys.argv) < 2:
    print 'Usage: %s [entry_name_1] [entry_name_2] ...'

entries = sys.argv[1:]
path = os.getcwd()

for lhs, rhs in combinations(entries, 2):
    raw_lhs = open(os.path.join(path, 
        '%s_model.txt' % lhs)).read().strip().splitlines(0)
    raw_rhs = open(os.path.join(path, 
        '%s_model.txt' % rhs)).read().strip().splitlines(0)

    ncoarse_lhs, atoms_lhs = find_meta(raw_lhs)
    ncoarse_rhs, atoms_rhs = find_meta(raw_rhs)

    m_source = """compareMotion('%(pdb_lhs)s', ...
    '%(pdb_rhs)s', ...
    %(ncoarse_lhs)s, %(ncoarse_rhs)s, ...
    [%(atoms_lhs)s], [%(atoms_rhs)s], ...
    7, 13)
""" % {
            'pdb_lhs': os.path.join(path, '%s_mixed_cg.pdb' % lhs),
            'pdb_rhs': os.path.join(path, '%s_mixed_cg.pdb' % rhs),
            'ncoarse_lhs': ncoarse_lhs,
            'ncoarse_rhs': ncoarse_rhs,
            'atoms_lhs': atoms_lhs,
            'atoms_rhs': atoms_rhs
        }

    fd = open("compare_%s_%s.m" % (lhs, rhs) , 'w+')
    fd.write(m_source)
    fd.close()

