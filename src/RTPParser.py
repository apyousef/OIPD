import re
from collections import OrderedDict

_RESIDUE_SECTION = re.compile(r'^\[\s+([A-Za-z0-9]*)\s+\]')
_ATOMS_SECTION = re.compile(r'^\s+\[\s+atoms\s+\]')
_ATOM_LINE = re.compile(r'[\s\t]+([A-Z0-9]+)[\s\t]+([A-Z0-9]+)[\s\t]+([-+]?[0-9]*\.?[0-9]+)[\s\t]+(\d+).*\n$')
_BONDS_SECTION = re.compile(r'^\s+\[\s+bonds\s+\]')
_BOND_LINE = re.compile(r'[\s\t]+([-+]?[A-Z0-9]+)[\s\t]+([-+]?[A-Z0-9]+)')
_GENERIC_SECTION = re.compile(r'^\s+\[\s+.*\s+\]')
FILEPATH = '/Users/andrewyousef/Documents/Ideas/ProteinModeling/OIPD/charmm27.ff/aminoacids.rtp'

class _AbstractParser(object):
    def __init__(self, _rtp):
        self._rtp = _rtp

    def parse(self, line):
        raise NotImplemented("Subclasses must implement instance method parse(self, line)")

    def next_parser(self, line):
        raise NotImplemented("Subclasses must implement instance method next_parser(self, line)")

class _InitialParser(_AbstractParser):
    def parse(self, line):
        pass

    def next_parser(self, line):
        if _RESIDUE_SECTION.match(line):
            return _ResidueParser(self._rtp)
        return self

class _ResidueParser(_AbstractParser):
    def parse(self, line):
        match = _RESIDUE_SECTION.match(line)
        if match:
            name = _RESIDUE_SECTION.match(line).groups()[0]
            residue = _Residue(name)
            self._rtp.add_residue(residue)

    def next_parser(self, line):
        if _ATOMS_SECTION.match(line):
            return _AtomParser(self._rtp)
        return self

class _AtomParser(_AbstractParser):
    def parse(self, line):
        match = _ATOM_LINE.match(line)
        if match:
            at_name = match.groups()[0]
            at_type = match.groups()[1]
            charge = float(match.groups()[2])
            charge_type = int(match.groups()[3])
            atom = _Atom(at_name, at_type, charge, charge_type)
            self._rtp.current_residue.add_atom(atom)

    def next_parser(self, line):
        if _BONDS_SECTION.match(line):
            return _BondParser(self._rtp)
        elif _RESIDUE_SECTION.match(line):
            return _ResidueParser(self._rtp)
        elif _GENERIC_SECTION.match(line):
            return _GenericParser(self._rtp)
        return self

class _BondParser(_AbstractParser):
    def parse(self, line):
        match = _BOND_LINE.match(line)
        if match:
            atom1 = match.groups()[0]
            atom2 = match.groups()[1]
            bond = (atom1, atom2)
            self._rtp.current_residue.add_bond(bond)

    def next_parser(self, line):
        if _RESIDUE_SECTION.match(line):
            return _ResidueParser(self._rtp)
        elif _GENERIC_SECTION.match(line):
            return _GenericParser(self._rtp)
        return self

# might get rid of this or _InitialParser since they seem to be identical
class _GenericParser(_AbstractParser):
    def parse(self, line):
        pass

    def next_parser(self, line):
        if _RESIDUE_SECTION.match(line):
            return _ResidueParser(self._rtp)
        return self

class _Atom(object):
    def __init__(self, at_name, at_type, charge, charge_group):
        self._at_name = at_name
        self._at_type = at_type
        self._charge = charge
        self._charge_group = charge_group

    @property
    def at_name(self):
        return self._at_name

    @property
    def at_type(self):
        return self._at_type

    @property
    def charge(self):
        return self._charge

# class _Bond(object):
#     def __init__(self, atom1, atom2):
#         self._atom1 = atom1
#         self._atom2 = atom2

class _Residue(object):
    def __init__(self, name):
        self._name = name
        self._atoms = dict()
        self._bonds = list()

    def add_atom(self, atom):
        self._atoms[atom.at_name] = atom

    def add_bond(self, bond):
        self._bonds.append(bond)

    @property
    def name(self):
        return self._name

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

class _RTP_internal(object):
    def __init__(self):
        self._residues = OrderedDict()

    def add_residue(self, residue):
        self._residues[residue.name] = residue

    @property
    def residues(self):
        return self._residues

    @property
    def current_residue(self):
        key = self._residues.keys()[-1]
        return self._residues.get(key)

class RTP(object):
    def __init__(self, filepath):
        self._rtp = _RTP_internal()
        current_parser = _InitialParser(self._rtp)
        with open(filepath, 'r') as f:
            for line in f:
                current_parser = current_parser.next_parser(line)
                current_parser.parse(line)

    @property
    def residues(self):
        return self._rtp._residues

r = RTP(FILEPATH)

print r.residues.get('ALA').atoms
print r.residues.get('ASPP').atoms['CG'].charge
print r.residues.get('ASPP').bonds