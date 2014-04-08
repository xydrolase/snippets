#!/usr/bin/env python

""" mixedcg.py
    A PyMOL extension script for building mixed coarse grained model. 

    Author: Xin Yin <xinyin at iastate dot edu>

    Biopython dependencies:
        This script requires Biopython module for communicating with NCBI
        server. You should first install Biopython onto your computer, and then
        copy the Bio directory installed at 
            /path/to/python/lib/site-packages/
        (or create a symbolic link of the Bio directory) to 
            /path/to/pymol/ext/lib/site-packages

        Try
            > import Bio
        in PyMOL to examine if you have correctly installed Biopython (for
        PyMOL).
        
    Example usage:
        
        In this example we are intended to prepare PDB files for MAVEN to run 
        mixed coarse-grained ANM, say the source PDB entry is 1A7S.

        > fetch 1a7s
        > run /path/to/mixedcg.py   # this step loads this script and registered
                                    # the cmd.mixed_cg() method

        > cmd.mixed_cg('1a7s',      # name of structure
                       around=10,   # use 10 angstroms as cutoff radius
                                    # for the PyMOL "object around radius"
                                    # selector
                       byres=True,  # Use by residue selecting mode
                       save=True    # Save the PDB files immediately
                       )

        After running this, you will have 4 new files, namely,
            1a7s_(fine|coarse)_grained.pdb
            1a7s_mixed_cg.pdb
            1a7s_model.txt
        
        You can then import 1a7s_mixed_cg.pdb to MAVEN for running mixed-cg ANM. 

        1a7s_model.txt contains the the number of atoms in each PDB file you
        sent to concatenate and also the atom indices for the active site
        residues that you will need for motion comparisons. 

        Example:

                    Atoms in each file:
                    1nes_coarse_grained.pdb - 180
                    1nes_fine_grained.pdb - 460

                    Atoms in active site residues:
                    263:272,353:357,476:480

        You will need all these numbers for running ANM and the dynamics
        comparison.

        --------

        For PDB entry 1NES, the annotated active sites from NCBI are incorrect.
        You can then explicitly specify the active sites by passing the
        argument `active_sites` to cmd.mixed_cg().

        Also, 1NES PDB file has 3 chains, namely chain E, I, J, which are the
        enzyme and two ligands respectively. The script will identify the
        longest chain as the protein for simulation. In cases you want to have 
        fine grained atoms to be selected around the ligands rather than the active
        sites, you can use additional parameter `ligand_chain_id`.

        > cmd.mixed_cg('1nes', ligand_chain_id='IJ', 
            active_sites=[57, 102, 195], save=True)

        Note that some parameters are omitted for the default values.

    Issues:
        I have not tried this on dimers or tetramers. Behavior for running this
        script on polymers might be odd.
"""

from Bio import Entrez

import os
import sys
import urllib
import re
 
# Identify yourself with the NCBI server
Entrez.email = "example@abc.com"

class EntrezRetriever:
    def __init__(self, db='gene'):
        self.db = db
        self.search_feature_header = re.compile(' {5}[A-za-z]+').search
        self.search_active_site_type = re.compile(r'/site_type="active').search
        self.search_active_sites = re.compile(r'order\((.+)\)').search

    def search_term(self, term, max_results=10):
        handle = Entrez.esearch(self.db, retmax=max_results, term=term)
        try:
            record = Entrez.read(handle)
        except RuntimeError, e:
            print "An error occurred during searching database."
            print "ERROR: %s" % e
            sys.exit(1)

        if record['Count'] > 0:
            return record['IdList']
        else:
            return []

    def retrieve_annotation(self, id_list):
        """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
        submit the data to NCBI) and esummary to retrieve the information. 
        Returns a list of dictionaries with the annotations."""
     
        request = Entrez.epost(self.db, id=",".join([str(id) for id in id_list]))
        try:
            result = Entrez.read(request)
        except RuntimeError, e:
            #FIXME: How generate NAs instead of causing an error with invalid IDs?
            print "An error occurred while retrieving the annotations."
            print "The error returned was %s" % e
            sys.exit(-1)
     
        webEnv = result["WebEnv"]
        queryKey = result["QueryKey"]
        data = Entrez.esummary(db=self.db, webenv=webEnv, query_key = queryKey)
        annotations = Entrez.read(data)
     
        return annotations

    def retrieve_features(self, id, feature):
        raw_metadata = self.strip_tags(urllib.urlopen(
            'http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?val=%(id)s&db=\
%(db)s&retmode=html' % {'id' : id, 'db': self.db}).read()).splitlines(True)

        lines_with_headers = [lineno for lineno in xrange(len(raw_metadata))
                if self.search_feature_header(raw_metadata[lineno])]
        lines_with_features = filter(
            lambda lineno: raw_metadata[lineno].startswith('     %s' % feature),
            lines_with_headers)

        features_selected = []
        for site_lineno in lines_with_features:
            index = lines_with_headers.index(site_lineno)
            if index < len(raw_metadata) - 1:
                features_selected.append(raw_metadata[lines_with_headers[index]:\
                        lines_with_headers[index+1]])
            elif index == len(raw_metadata):
                features_selected(raw_metadata[lines_with_features[index]:])

        return features_selected

    def find_active_sites(self, all_sites, offset):
        # Only return the first one?
        for site in all_sites:
            if self.search_active_site_type(site[1]):
                return map(lambda acs:
                        int(acs)+offset,
                        self.search_active_sites(site[0]).group(1).split(','))

        return None

    def strip_tags(self, html):
        return re.sub('</?[^<>]+?>', '', html)

def construct_mixed_cg(object_id, around_with, selector_peptide, radius, byres):
    cmd.select('%s_fine_grained' % object_id,
            '%(byres)s %(object)s_%(with)s \
around %(radius)d and %(object)s and %(chain)s and not resn HOH+SO4' % {
                'object': object_id,
                'chain': selector_peptide,
                'with': around_with,
                'radius': radius,
                'byres': ['', 'byres'][byres]
            })

    cmd.select('%s_coarse_grained' % object_id,
            '%(object)s and %(chain)s and not (%(object)s_%(with)s \
or %(object)s_fine_grained) and name CA' % {
                'object': object_id,
                'with': around_with,
                'chain': selector_peptide,
            })

    # Now try to re-visualize the molecule
    cmd.hide('everything', object_id)
    for color, shape, parts in [('red', 'sticks', around_with),
            ('marine', 'sticks', 'fine_grained'),
            ('gray', 'spheres', 'coarse_grained')]:

        selector = '%(object)s_%(parts)s' % {
                'object': object_id,
                'parts': parts
            }

        cmd.color(color, selector)
        cmd.show(shape, selector) 

def save_model(object_id, around_with, path=''):
    """Save the mixed cg model to PDB files at given path."""
    for suffix in ['coarse_grained', 'fine_grained']:
        pdb = os.path.join(path,
                '%(object)s_%(suffix)s.pdb' % {
                    'object': object_id,
                    'suffix': suffix
                })

        if suffix == 'fine_grained':
            selector = '%(object)s_fine_grained or %(object)s_%(with)s' % {
                    'object': object_id,
                    'with': around_with
                    }
        else:
            selector = '%s_coarse_grained' % object_id

        cmd.save(pdb, selector)

def concatenate_pdbs(pdbs, active_sites, path=''):
    pdbs = [os.path.join(path, 
        [''.join((f, '.pdb')), f][f.endswith('.pdb')]) 
        for f in pdbs]

    object_id = [f[:f.find('_')] for f in pdbs][0]

    atom_serial = 1
    prev_serial = 1

    active_site_atoms = []
    active_sites = map(str, active_sites)    # for string comparison

    atom_in_files = []
    output = []

    for coords_file in pdbs:
        if os.path.exists(coords_file):
            records = open(coords_file).read().splitlines(False)
            
            for rec in records:
                if rec[:6] in ('ATOM  ', 'HETATM'):
                    # Split coordinate record into non-trivial data and spaces
                    # so we can deal with specific columns.
                    coord_columns = re.split(r'\s+', rec)
                    coord_spaces = re.findall(r'\s+', rec)

                    # We will now reset all atom indices to sequential order.
                    # To maintain the alignment of the coordinate records, we
                    # need to compensate the alternations in atom indices by
                    # adjusting spaces
                    space_offset = len(coord_columns[1]) - len(str(atom_serial))
                    coord_spaces[0] = ' ' * (len(coord_spaces[0]) + space_offset)
                    # Reset the serial
                    coord_columns[1] = str(atom_serial)

                    # We need to record all the (modified) atom indices that
                    # corresponds to active site residues, which will be
                    # compared later.
                    if coord_columns[5] in active_sites:
                        active_site_atoms.append(atom_serial)

                    atom_serial += 1

                    # Regroup all data back into string
                    coord_all = []
                    for col, space in zip(coord_columns, coord_spaces):
                        coord_all.extend([col, space])

                    output.append(''.join(coord_all))

            atom_in_files.append('%s - %d' % (
                coords_file, atom_serial - prev_serial
                ))

            prev_serial = atom_serial
        
        fd = open(os.path.join(path, '%s_mixed_cg.pdb' % object_id), 'w+')
        fd.write('\n'.join(output))
        fd.close()

        # Try to represent the indices for atoms in MATLAB compatible vector
        if active_site_atoms:
            _atoms = []
            atom_count = len(active_site_atoms) - 1
            reduce(lambda x, y: _atoms.append(
                [str(active_site_atoms[x]), '-'][y < atom_count and\
                        active_site_atoms[y]-active_site_atoms[x]==1]) or y, 
                xrange(len(active_site_atoms)))

            vec_atoms = re.sub('(-+)(\d+)', 
                    lambda x: "%d:%s," % (
                        int(x.group(2))-len(x.group(1)),
                        x.group(2)),
                    ''.join(_atoms))[:-1]
        else:
            vec_atoms = ''

        fd = open(os.path.join(path, '%s_model.txt' % object_id), 'w+')
        fd.write("Atoms in each file:\n" + '\n'.join(atom_in_files) +\
                '\n\nAtoms in active site residues:\n')
        fd.write(vec_atoms + '\n')
        fd.close()

def mixed_cg(object_id, ligand_chain_id=None, active_sites=None,
        around=10, offset=0, byres=False, save=False):
    entrezer = EntrezRetriever(db='protein')

    chains = sorted([(len(cmd.get_model('%(object)s and chain %(chain)s' %
        {'object': object_id, 'chain': chain_id}).atom), chain_id) 
        for chain_id in cmd.get_chains(object_id)],
        lambda x, y: cmp(y[0], x[0]))

    selector_peptide = 'chain %s' % chains[0][1]
    entrez_protein_id = '%(object)s_%(chain)s' % {
            'object': object_id, 'chain': chains[0][1]}

    if not active_sites:
        print 'Searching NCBI using', entrez_protein_id

        id_list = entrezer.search_term(entrez_protein_id)
        if id_list:
            # Grab all annotations metadata categorized as functional sites
            print 'Found protein entry', id_list[0]
            active_sites = entrezer.find_active_sites(
                    entrezer.retrieve_features(id_list[0], 'Site'), offset)
            if active_sites:
                print 'Found active sites', \
                        ','.join(map(str, active_sites)), 'for', object_id

    if not ligand_chain_id and active_sites:
        # Active sites found, try build mixed cg model by selecting atoms
        # close to proximity of active sites as fine grained atoms.
        cmd.select('%s_active_sites' % object_id,
                '(%(object)s and %(chain)s) and resi %(residues)s' % {
                    'object': object_id, 
                    'chain': selector_peptide,
                    'residues': '+'.join(map(str, active_sites))
                })

        construct_mixed_cg(object_id, 'active_sites',
                selector_peptide, around, byres)
    else:
        # User specified chain(s) as ligand(s), build mixed cg model based on
        # that.
        selector_ligands = ' or '.join(map(lambda arg: 'chain %s' % arg,
            list(ligand_chain_id)))
        ligands = cmd.select('%s_ligands' % object_id,
                '%(object)s and (%(ligand_chains)s) and not resn HOH+SO4' % {
                    'object': object_id,
                    'ligand_chains': selector_ligands
                })

        construct_mixed_cg(object_id, 'ligands', selector_peptide, around,
                byres)

    if save:
        # Save the PDB files to current directory, and perform concatenation
        save_model(object_id, 
                ['ligands', 'active_sites'][ligand_chain_id == None])
        concatenate_pdbs(('%s_coarse_grained.pdb' % object_id, 
                '%s_fine_grained.pdb' % object_id),
                active_sites)

        print 'Mixed CG model PDB files saved to your current directory.'
        
# Insert the function into cmd namespace
# NOTE: We avoid using cmd.extend() functionality here as suggested by PyMOL
# scripting tutorial. Calling cmd.mixed_cg() brings more flexibility in using
# optional parameters and consistency in avoiding weird errors. (As PyMOL
# treats every argument, even those supposed to be integers, as strings)

cmd.mixed_cg = mixed_cg
