#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
protstruct.py
Utility module that defines classes related to protein structures

Author: Xin Yin <xinyin at iastate dot edu>

DESCRIPTION
This module implements interface to protein structures as well as
consisting atoms/residues.
"""

__version__ = 0.1
__author__ = "Xin Yin"

from copy import copy

import math
import re
import os

import numpy as np

import geometry

sign = lambda x: -1 if x < 0 else 1

class Atom:
    name = ''
    serial = 0 
    altLoc = ''
    resSeq = 0
    iCode = ''
    resName = ''
    chainID = 0
    x = 0.0
    y = 0.0
    z = 0.0
    occupancy = 0.0
    tempFactor = 0.0
    element = '' 
    charge = '' 

    # As per PDB format 3.3
    atom_rec_columns = [
            #(slice(0, 6), 'recName', str),
            (slice(6, 11), 'serial', int),
            (slice(12, 16), 'name', str),
            (slice(16, 17), 'altLoc', str),
            (slice(17, 20), 'resName', str),
            (slice(21, 22), 'chainID', str),
            (slice(22, 26), 'resSeq', int),
            (slice(26, 27), 'iCode', str),
            (slice(30, 38), 'x', float),
            (slice(38, 46), 'y', float),
            (slice(46, 54), 'z', float),
            (slice(54, 60), 'occupancy', float),
            (slice(60, 66), 'tempFactor', float),
            (slice(76, 78), 'element', str),
            (slice(78, 80), 'charge', str),
            ]

    def __init__(self, dummy, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def identify(self):
        return '{0:5d} {1:s} ({2:s} {3:4d})'.format(
                self.serial, self.name.rjust(4), self.resName, 
                self.resSeq)

    def __repr__(self):
        return '{0:5d} {1:s} ({2:s}) <{3:0.4f}, {4:0.4f}, {5:0.4f}>'.format(
                self.serial, self.name.rjust(4), self.resName,
                self.x, self.y, self.z)

    def __str__(self):
        return 'ATOM  {0:5d} {1}{2}{3} {4}{5:4d}{6}   {7:8.3f}{8:8.3f}\
{9:8.3f}{10:6.2f}{11:6.2f}          {12} {13}\n'.format(
         self.serial, self.name.ljust(4), self.altLoc.ljust(1), self.resName,
         self.chainID, self.resSeq, self.iCode.ljust(1),
         self.x, self.y, self.z, self.occupancy, self.tempFactor,
         self.element.rjust(2), self.charge.rjust(2))

    def __sub__(self, atom_b):
        """Compute the vector corresponding to (atom_a - atom_b) in
        Cartesian space, where `atom_a` here is the `self` object.

        This operation returns a `numpy.array` type."""

        return np.array([self.x - atom_b.x, self.y - atom_b.y,
            self.z - atom_b.z])

    def direction_to(self, atom_b):
        """Returns a unit vector corresponding to the direction from
        atom_a to atom_b, where `atom_a` here is the `self` object."""

        _dir = np.array([atom_b.x - self.x, atom_b.y - self.y,
            atom_b.z - self.z])

        return _dir / np.linalg.norm(_dir)

    @classmethod
    def distance(self, atom_a, atom_b):
        """Compute distance between any two atoms. If the two atoms are
        adjacent on the backbone of a protein, this distance coincides with
        the bond length."""

        return np.sqrt(np.sum((atom_b - atom_a) ** 2))

    @classmethod
    def angle(self, atom_a, atom_b, atom_c):
        """Compute bond angle between three adjacent atoms."""

        v_ba = atom_a - atom_b
        v_bc = atom_c - atom_b
        return math.acos(np.vdot(v_ba, v_bc) / (
            np.linalg.norm(v_ba) * np.linalg.norm(v_bc)))

    @classmethod
    def dihedral(self, atom_a, atom_b, atom_c, atom_d):
        """Compute the torsional or dihedral angle between two planes
        <a-b-c> and <b-c-d>. The sign of the angle is determined by how
        the bond <a-b> is rotated relative to the normal of <b-c-d> plane,
        using right-handness rule."""

        # normal of two planes
        nv_abc = np.cross(atom_a - atom_b, atom_c - atom_b)
        nv_bcd = np.cross(atom_b - atom_c, atom_d - atom_c)

        # the angle between two normals
        normal_angle = math.acos(np.vdot(nv_abc, nv_bcd) / (
            np.linalg.norm(nv_abc) * np.linalg.norm(nv_bcd)))

        # -- cis or trans:
        #    compute the cross product between (b-a) and the normal of 
        #    plane <b-c-d>: if positive: handness is correct, assign 
        #    positive sign to the dihedral angle, otherwise negative.
        normal_angle *= sign(np.vdot(atom_a - atom_b, nv_bcd))

        return normal_angle# / math.pi * 180

    @classmethod
    def parse(self, pdb_str):
        """Convert a string literal ATOM record into an instance of class
        Atom."""
        atom_desp = ((k, conv(pdb_str[_slc].strip()))
                for _slc, k, conv in self.atom_rec_columns)
        return Atom(None, **dict(atom_desp))

    def dist_to(self, atom_b):
        return Atom.distance(self, atom_b)

    def xyz(self):
        return np.array([self.x, self.y, self.z])

    def relocate(self, xyz):
        self.x, self.y, self.z = list(xyz)

class Residue:
    name = None
    sequence = 0
    atom_index = {}
    atoms = []

    def __init__(self, name, seq):
        self.name = name
        self.sequence = seq

    def __repr__(self):
        return '{0} {1:4d}'.format(self.name.ljust(4), self.sequence)

    def __str__(self):
        return self.__repr__()
    
    def __setattr__(self, key, value):
        """Wraps around magic method __setattr__ to define a trigger such that
        each time user assigns a new atom list, call a callback function to
        update atom indexing."""
        self.__dict__[key] = value
        if key == 'atoms':
            self.index_atoms()

    def __getitem__(self, key):
        if type(key) is int:
            return self.atoms[key]
        elif type(key) is str:
            return self.atom_index[key]

    def index_atoms(self):
        self.atom_index = {}
        for atom in self.atoms:
            self.atom_index[atom.name] = atom


class Protein:
    """A simple PDB file parser."""

    pdb_file = None
    pdb_content = None
    atoms = []
    backbone = []
    residues = []

    backbone_filter = set(['C', 'CA', 'N'])

    # http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
    sidechain_di_atoms = {
        'CYS': [
            ['N','CA','CB','SG']
        ], 
        'ASP': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','OD1']
        ],
        'SER': [
            ['N','CA','CB','OG']
        ], 
        'GLN': [
            ['N','CA','CB','CG'],
            ['CA','CB','CG','CD'],
            ['CB','CG','CD','OE1']
        ],
        'LYS': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','CD'],
            ['CB','CG','CD','CE'], 
            ['CG','CD','CE','NZ']
        ], 
        'ILE': [
            ['N','CA','CB','CG1'], 
            ['CA','CB','CG1','CD1']
        ], 
        'PRO': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','CD']
        ], 
        'THR': [
            ['N','CA','CB','OG1']
        ], 
        'PHE': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','CD1']
        ],
        'ASN': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','OD1']
        ], 
        'MET': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','SD'],
            ['CB','CG','SD','CE']
        ], 
        'HIS': [
            ['N','CA','CB','CG'],
            ['CA','CB','CG','ND1']
        ],
        'LEU': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','CD1']
        ], 
        'ARG': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','CD'], 
            ['CB','CG','CD','NE'], 
            ['CG','CD','NE','CZ'],
            ['CD','NE','CZ','NH1']
        ], 
        'TRP': [
            ['N','CA','CB','CG'],
            ['CA','CB','CG','CD1']
        ],
        'VAL': [
            ['N','CA','CB','CG1']
        ], 'GLU': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','CD'], 
            ['CB','CG','CD','OE1']
        ],
        'TYR': [
            ['N','CA','CB','CG'], 
            ['CA','CB','CG','CD1']
        ]}

    def __init__(self, pdb_file):
        if os.path.exists(pdb_file):
            self.pdb_file = pdb_file
            with open(pdb_file) as f:
                self.pdb_content = f.readlines()

            self._extract_atoms()
            self._extract_residues()
            self._compute_sidechains()

    def __len__(self):
        """Return the number of atoms in the protein."""
        return len(self.atoms)

    def __getitem__(self, index):
        if type(index) is int:
            return self.atoms[index]
        else:
            raise TypeError("Atom index has to be an integer.")

    def _compute_sidechains(self):
        """Identifying covalent bonds using empirical rules used by RasMol.
        See http://www.umass.edu/microbio/rasmol/rasbonds.htm for discussion
        on the such rules.

        After identifying bonds for the atom pairs within each residue, the
        sidechains are identified as the set of atoms that are (directly
        or indirectly through other bonds) connected to the backbone atoms
        (N, CA and C).

        This function prepares sidechain information so that when a bond is
        rotated for a given torsional angle, the sidechain atoms can be
        readily identified and rotated as well.
        """

        def rel_internal_coords(sc, a):
            """Compute internal coordinates of a side chain atom `sc` 
            relative to the backbone atom `a`."""
            if a.backbone_index >= len(self.backbone) - 2:
                return (0, 0, 0)

            b, c = self.backbone[a.backbone_index+1:a.backbone_index+3]
            return (Atom.distance(sc, a), Atom.angle(sc, a, b),
                    Atom.dihedral(sc, a, b, c))

        for res in self.residues:
            atoms_done = set([atom for atom in res.atoms \
                    if atom.name in self.backbone_filter])

            # atoms to build bonds
            atoms_todo = set(res.atoms) - atoms_done

            while atoms_todo:
                connected_atoms = set()
                for a in atoms_todo:
                    for b in atoms_done:
                        # --- empirical rules used by RasMol ---
                        if a.element == 'H' or b.element == 'H':
                            if 0.4 < Atom.distance(a, b) < 1.2:
                                connected_atoms.add(a)
                                a.connected_to = b.connected_to
                                a.connected_to.sidechain.append(a)
                                a.int_coords = rel_internal_coords(a,
                                        a.connected_to)
                        else:
                            if 0.4 < Atom.distance(a, b) < 1.9:
                                connected_atoms.add(a)
                                a.connected_to = b.connected_to
                                a.connected_to.sidechain.append(a)
                                a.int_coords = rel_internal_coords(a,
                                        a.connected_to)

                atoms_todo -= connected_atoms
                atoms_done |= connected_atoms

    def _extract_residues(self):
        res_atoms = []
        res_seq = self.atoms[0].resSeq

        for atom in self.atoms:
            if atom.resSeq != res_seq:
                new_res = Residue(res_atoms[0].resName, res_seq)
                new_res.atoms = copy(res_atoms)
                self.residues.append(new_res)

                res_atoms = []

            res_seq = atom.resSeq
            res_atoms.append(atom)

        new_res = Residue(res_atoms[0].resName, res_seq)
        new_res.atoms = copy(res_atoms)
        self.residues.append(new_res)

    def _extract_atoms(self):
        """Refer to:
            http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
        for PDB format specifications."""

        regex_atom = re.compile("^ATOM")
        atom_entries = [line for line in self.pdb_content \
                if regex_atom.search(line)]

        self.atoms = [Atom.parse(rec) for rec in atom_entries]
        self.backbone = [atom for atom in self.atoms \
                if atom.name in self.backbone_filter]

        for index, atom in enumerate(self.backbone):
            atom.backbone_index = index
            atom.connected_to = atom
            atom.sidechain = []

    def write_pdb(self, file):
        with open(file, 'w') as f:
            atom_idx = [idx for idx in range(len(self.pdb_content)) \
                    if self.pdb_content[idx].startswith('ATOM')]
            f.writelines(self.pdb_content[:min(atom_idx)])
            f.writelines(str(atom) for atom in self.atoms)
            f.writelines(self.pdb_content[max(atom_idx)+1:])

    def compute_sidechain_dihedrals(self):
        dihedrals = []
        for res in self.residues:
            di_def_atoms = self.sidechain_di_atoms.get(res.name, None)
            if di_def_atoms is None:
                continue

            res_dihedrals = []
            for angle_atoms in di_def_atoms:
                atom_list = [res[name] for name in angle_atoms]
                res_dihedrals.append(Atom.dihedral(*atom_list))

            dihedrals.append((res, res_dihedrals))

        return dihedrals

    def cartesian_coords(self):
        """Return a Nx3 numpy.matrix within which each row corresponds to
        an atom in the protein backbone. The three columns represent the 
        x, y and z coordinates in space respectively."""

        return np.matrix([atom.xyz() for atom in self.backbone])

    def internal_coords(self):
        """Compute the internal coordinates for the polymer structure
        (backbone only), specified by 3N-1 bond lengths, 3N-2 bond angles
        and 3N-3 dihedral angles where N is the number of residues."""
        n_atoms = len(self.backbone)

        bond_lens, bond_angles, dihedrals = zip(*[
            (
                Atom.distance(*self.backbone[p:p+2]),
                Atom.angle(*self.backbone[p:p+3]),
                Atom.dihedral(*self.backbone[p:p+4])
            ) for p in xrange(0, n_atoms - 3)])

        bond_angles = list(bond_angles) + [
                Atom.angle(*self.backbone[n_atoms-3:])
                ]
        bond_lens = list(bond_lens) + [
                Atom.distance(*self.backbone[n_atoms-3:n_atoms-1]),
                Atom.distance(*self.backbone[n_atoms-2:n_atoms])
                ]

        return bond_lens, bond_angles, dihedrals

    def rearrange_internal(self, bond_lens, bond_angles, dihedrals):
        """Rearrange the structe of the protein using the specified
        internal coordinates of the atoms. All Cartesian coordinates 
        will be recomputed based on the internal coordinates.

        To resolve identifiability issues in Cartesian space, this
        function will use the Cartesian coordinates (x, y and z) of
        the last three atoms in the protein to build a reference plane.
        Hence the extra 9 degrees of freedom of the external coordinates
        can be uniquely specified."""

        _, __, di_orig = self.internal_coords()

        n_bb_atoms = len(self.backbone)

        if len(bond_lens) != n_bb_atoms - 1 or \
                len(bond_angles) != n_bb_atoms - 2 or \
                len(dihedrals) != n_bb_atoms - 3:
            raise ValueError("Invalid internal coordinates.")

        # Compute the cartesian coordinates of atom A given the adjacent
        # backbone atoms B, C and D at its C-terminus. 
        for atom_index in range(n_bb_atoms-4, -1, -1):
            b, c, d = self.backbone[atom_index+1:atom_index+4]
            a = self.backbone[atom_index]

            dir_bc = b.direction_to(c)
            # the normal of plane B-C-D 
            normal_bcd = geometry.normal_p3(b, c, d)

            tor_angle = dihedrals[atom_index]
            # vector from B to A, using Rodrigues formula to compute
            # vector rotations related to bond angles and dihedral angles.
            v_ba_cis = geometry.rotate_about(dir_bc, -normal_bcd,
                    bond_angles[atom_index]) * bond_lens[atom_index]
            v_ba = geometry.rotate_about(v_ba_cis, -dir_bc, tor_angle)

            a.relocate(b.xyz() + v_ba)

            # now compute the coordinates for sidechain atoms based on 
            # its relative position to the backbone atom.
            normal_abc = geometry.normal_p3(a, b, c)
            dir_ab = a.direction_to(b)

            for sc_atom in a.sidechain:
                sc_bl, sc_ba, sc_di = sc_atom.int_coords
                v_asc_cis = geometry.rotate_about(dir_ab, -normal_abc,
                        sc_ba) * sc_bl
                v_asc = geometry.rotate_about(v_asc_cis, -dir_ab, sc_di) 

                sc_atom.relocate(a.xyz() + v_asc)

    def translate(self, x, y, z):
        pass

    def rotate(self, phi, psi, omega):
        pass
