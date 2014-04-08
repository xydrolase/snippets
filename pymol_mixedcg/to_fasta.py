#!/usr/bin/env python

import re

AA_map = {
    "ALA":"A",
    "ARG":"R",
    "ASN":"N",
    "ASP":"D",
    "CYS":"C",
    "GLU":"E",
    "GLN":"Q",
    "GLY":"G",
    "HIS":"H",
    "ILE":"I",
    "LEU":"L",
    "LYS":"K",
    "MET":"M",
    "PHE":"F",
    "PRO":"P",
    "SER":"S",
    "THR":"T",
    "TRP":"W",
    "TYR":"Y",
    "VAL":"V",
}

def to_fasta(objs, sites):
    if isinstance(objs, str):
        objs = [objs]

    for entry in objs:
        cmd.select("%s_sites" % entry, "%(object)s and resi %(sites)s" % {
            'object': entry, 
            'sites': '+'.join(map(str, sites))
        })

        cmd.select("%s_close_residues" % entry,
                   "byres %(object)s_sites around 7.5 and %(object)s" % {
                       'object': entry
                })

        selections = "(%(object)s_close_residues or %(object)s_sites)\
and name CA" % {
                     'object': entry
                }

        cmd.save("%s_res.pdb" % entry, selections)

        records = map(lambda y: re.split('\s+', y[17:]), filter(
            lambda x: x.startswith('ATOM'),
            open("%s_res.pdb" % entry, "r").read().splitlines(0)
        ))

        aa_chars = map(lambda aa: AA_map[aa[-3:]],
                       [res
                        for res, index in 
                            [(e[0], e[2]) for e in records]
                        ])

        fd = open("%s_res.fasta" % entry, 'w+')
        fd.write('''\
>%s
%s
''' 
                 % (entry, ''.join(aa_chars)))
        fd.close()


cmd.sites2fasta = to_fasta


