import numpy as np


class FakeParent:
    def __init__(self, id):
        self.id = str(id)

    def __repr__(self):
        return f'\nFakeParent id={self.id}'


class FakeResidue():
    def __init__(self, dict_atoms=None, name=None, _id=()):
        self.dict_atoms = dict_atoms
        self.name = name
        self.id = ['resname', _id, ' ']
        self.resname = name
        self.dict_data = {}


    def __getitem__(self, key):
        return self.dict_atoms[key]

    def get_atoms(self):
        return self.dict_atoms.values()

    def get_resname(self):
        return self.name

    @property
    def xyz(self):
        return 0,0,0


    @staticmethod
    def get_parent():
        return FakeParent('X')


class FakeChain:
    def __init__(self, id):
        self.id = id

    def __repr__(self):
        return f'<Chain id={self.id}>'


class FakeChild:
    def __init__(self, list_chains, id=None):
        if id is None:
            self.id = '<Model id=0>'
        else:
            self.id = f'<Model id={id}>'
        child_list = []
        for chain in list_chains:
            child_list.append(FakeChain(chain))
        self.child_list = child_list

    def __repr__(self):
        return f'{self.id}'

    def detach_child(self, chain_ID):
        for child in self.child_list:
            if child.id == chain_ID:
                self.child_list.remove(child)


class FakePDB:

    def __init__(self, list_atoms, id=None, dict_pigments=None, set_atomic_charges=False, child_list=None,
                 list_chains=None):
        self.list_atoms = list_atoms
        self.id = id
        self.dict_pigments = dict_pigments
        if dict_pigments is None:
            self.dict_pigments = {}
        if list_chains is None:
            list_chains = ['X', 'Y', 'Z']
        if child_list is None:
            child_list = [FakeChild(list_chains=list_chains)]
        self.child_list = child_list

        if set_atomic_charges:
            self.set_atomic_charges(pdb)

    def get_models(self):
        return self

    def get_chains(self):
        return self.child_list

    def get_atoms(self):
        return self.list_atoms

    def transform(self, rot=np.eye(3), tran=np.zeros(3)):
        for atom in self.list_atoms:
            atom.transform(rot, tran)

    def get_residues(self):
        return self.child_list

    def __repr__(self):
        return f'\nFakePDB {self.id} (list_atoms = {self.list_atoms}, dict_pigments = {self.dict_pigments}'


class FakeAtom:

    def __init__(self, xyz, mass=1, id=0, pqr_charge=None, element='X'):
        self.coord = np.array(xyz, dtype=float)
        self.mass = mass
        self.id = id
        self.element = element
        self.pqr_charge = pqr_charge

    def __repr__(self):
        return f'\nFakeAtom {self.id} (element = {self.element}, coord = {self.coord}, mass = {self.mass}'

    def transform(self, rot=np.eye(3), tran=np.zeros(3)):
        self.coord = np.dot(self.coord, rot) + tran
