#%%
import sys
import math
import re
import os
from io import BytesIO, StringIO
import requests, collections
import _cobio as cb

#%%
class PdbAtom:
    def __init__(self, name=None, num=None, coords=None, is_het=None):
        self.name = "" if name is None else name
        self.num = -1 if num is None else num
        self.coords = [0, 0, 0] if coords is None else coords
        self.is_het = False if is_het is None else is_het

    def __getitem__(self, i):
        return self.coords[i]

    def __setitem__(self, i, v):
        self.coords[i] = v

    def __iter__(self):
        return iter(self.coords)

    @classmethod
    def from_dict(cls, dt):
        return PdbAtom(**dt)


class PdbResidue:
    def __init__(self, name=None, num=None, atoms=None):
        self.name = "" if name is None else name
        self.num = -1 if num is None else num
        self.atoms = [] if atoms is None else atoms

    def __iter__(self):
        return iter(self.atoms)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.atoms[key]
        else:
            for atom in self.atoms:
                if atom.name == key:
                    return atom
            print("Couldn't find atom '%s'" % key)
            write_residue(self)
            quit()

    @classmethod
    def from_dict(cls, dt):
        if "atoms" in dt:
            dt["atoms"] = [PdbAtom.from_dict(dt) for dt in dt["atoms"]]
        return PdbResidue(**dt)


class PdbChain:
    def __init__(self, name=None, residues=None):
        self.name = "" if name is None else name
        self.residues = [] if residues is None else residues

    def __iter__(self):
        return iter(self.residues)

    def __getitem__(self, key):
        if isinstance(key, int):
            return super(PdbChain, self).__getitem__(key)
        else:
            for residue in self:
                if residue.num == key:
                    return residue
            print("Couldn't find residue %s in chain %s" % (key, self.name))
            quit()

    @property
    def seq(self):
        return "".join(res.name.strip()[0] for res in self.residues)

    @property
    def atoms(self):
        return [atom for residue in self.residues for atom in residue.atoms]

    @classmethod
    def from_dict(cls, dt):
        if "residues" in dt:
            dt["residues"] = [PdbResidue.from_dict(dt) for dt in dt["residues"]]
        return PdbChain(**dt)


class PdbModel:
    def __init__(self, num=None, chains=None):
        self.num = -1 if num is None else num
        self.chains = [] if chains is None else chains

    def __iter__(self):
        return iter(self.chains)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.chains[key]
        else:
            for chain in self:
                if chain.name == key:
                    return chain
            print("Couldn't find chain %" % key)
            quit()

    @property
    def seq(self):
        return "".join(res.name.strip()[0] for chain in self.chains for res in chain)

    @property
    def residues(self):
        return [residue for chain in self for residue in chain]

    @property
    def atoms(self):
        return [atom for chain in self for residue in chain for atom in residue]

    @classmethod
    def from_dict(cls, dt):
        if "chains" in dt:
            dt["chains"] = [PdbChain.from_dict(dt) for dt in dt["chains"]]
        return PdbModel(**dt)


class Pdb:
    def __init__(self, name=None, models=None):
        self.name = "" if name is None else name
        self.models = [] if models is None else models

    def __iter__(self):
        return iter(self.models)

    def __getitem__(self, key):
        return self.models[key]

    @classmethod
    def from_dict(cls, dt):
        if "models" in dt:
            dt["models"] = [PdbModel.from_dict(dt) for dt in dt["models"]]
        return Pdb(**dt)


#%%
def read_pdb(filename="", content="", name=""):
    if filename != "":
        if name == "":
            name = os.path.splitext(os.path.basename(filename))[0]
        return PdbParser(open(filename), name).pdb
    if content != "":
        return PdbParser(StringIO(content), name).pdb


def read_cif_as_pdb(filename):
    return Pdb.from_dict(cb.read_cif_as_pdb(filename))


class ParsedLine:
    def __init__(self, line):
        self.atom_name = line[12:16].strip()
        self.res_name = line[17:20].strip()
        self.chain_name = line[20:22].strip()
        self.atom_num = int(line[6:11].strip())
        self.res_num = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        self.is_het = line.startswith("HETATM")


class PdbParser:
    def __init__(self, IO, name):
        self.atoms = []
        self.residues = []
        self.chains = []
        self.models = []
        self.oline = ""
        self.model_num = 1
        i = 0
        for line in IO:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                line = ParsedLine(line)
                if i > 0:
                    if line.chain_name != self.oline.chain_name:
                        self.add_chain()
                    elif line.res_num != self.oline.res_num or line.res_name != self.oline.res_name:
                        self.add_residue()
                self.atoms.append(PdbAtom(line.atom_name, line.atom_num, [line.x, line.y, line.z]))
                self.atoms[-1].is_het = line.is_het
                self.oline = line
                i += 1
            elif line[0:5] == "MODEL":
                self.add_model()
                g = re.split("\s+", line)
                if len(g) >= 2:
                    try:
                        self.model_num = int(g[1])
                    except ValueError:
                        pass
            elif line[0:6] == "ENDMDL":
                self.add_model()
            elif line[0:3] == "TER":
                self.add_chain()
        self.add_model()
        self.pdb = Pdb(name, self.models)

    def add_residue(self):
        if len(self.atoms) > 0:
            self.residues.append(PdbResidue(self.oline.res_name, self.oline.res_num, self.atoms))
            self.atoms = []

    def add_chain(self):
        self.add_residue()
        if len(self.residues) > 0:
            self.chains.append(PdbChain(self.oline.chain_name, self.residues))
            self.residues = []

    def add_model(self):
        self.add_chain()
        if len(self.chains) > 0:
            self.models.append(PdbModel(self.model_num, self.chains))
            self.chains = []
            self.model_num += 1


def print_object(o):
    for attr in dir(o):
        if not attr.startswith("__"):
            value = getattr(o, attr)
            print(f"{attr}: {value}")


def read_traj(filename):
    atom, residue, chain, model = PdbAtom(), PdbResidue(), PdbChain(), PdbModel()

    def add_residue(l):
        if len(residue.atoms) > 0:
            residue.name = l.res_name
            residue.num = l.res_num
            chain.residues.append(residue)
            return PdbResidue()
        else:
            return residue

    def add_chain(l):
        if len(chain.residues) > 0:
            chain.name = l.chain_name
            model.chains.append(chain)
            return PdbChain()
        else:
            return chain

    l, l_old = None, None
    for line in open(filename):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            l = ParsedLine(line)

            if len(residue.atoms) > 0:
                # if residue changes
                if l.res_num != l_old.res_num or l.res_name != l_old.res_name or l.chain_name != l_old.chain_name:
                    residue = add_residue(l_old)
                    # print(len(residue.atoms))
                # if chain changes
                if l.chain_name != l_old.chain_name:
                    chain = add_chain(l_old)

            # if model changes
            atom.num = l.atom_num
            atom.name = l.atom_name
            atom.coords = [l.x, l.y, l.z]
            residue.atoms.append(atom)
            atom = PdbAtom()

            l_old = l

        elif line.startswith("TER"):
            residue = add_residue(l_old)
            chain = add_chain(l_old)
        elif line.startswith("MODEL") or line.startswith("ENDMDL"):
            residue = add_residue(l_old)
            chain = add_chain(l_old)
            if len(model.chains) > 0:
                yield model
                model = PdbModel()
        elif line.startswith("END"):
            residue = add_residue(l_old)
            chain = add_chain(l_old)
            if len(model.chains) > 0:
                yield model
                return
    residue = add_residue(l_old)
    chain = add_chain(l_old)
    if len(model.chains) > 0:
        yield model
    return


#%%
def write_pdb(pdb, f=sys.stdout):
    model_index = 0
    for model in pdb:
        f.write("MODEL%6d\n" % (model_index + 1,))
        write_model(model, f)
        f.write("ENDMDL\n")
        model_index += 1
    f.write("END\n")


def pdb_line(
    record_name,
    atom_num,
    atom_name,
    residue_name,
    chain_name,
    residue_num,
    x,
    y,
    z,
    occupancy=1.0,
    temperature_factor=0.0,
    element="",
    charge="",
):
    if len(atom_name) == 2:
        atom_name += " "
    elif len(atom_name) == 1:
        atom_name += "  "
    return "%-6.6s%5d %4.4s %3.3s%2.2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12.12s%2s" % (
        record_name,
        atom_num,
        atom_name,
        residue_name,
        chain_name,
        residue_num,
        x,
        y,
        z,
        occupancy,
        temperature_factor,
        element,
        charge,
    )


def write_model(model, f=sys.stdout):
    atom_index = 0
    for chain in model:
        for residue in chain:
            for atom in residue:
                f.write(pdb_line("ATOM", atom_index + 1, atom.name, residue.name, chain.name, residue.num, atom[0], atom[1], atom[2]))
                atom_index += 1


def write_residue(residue):
    num_atom = 1
    for atom in residue:
        print(
            "ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n"
            % (num_atom, atom.name, residue.name, "A", residue.num, atom[0], atom[1], atom[2], 1.00, 0.00, atom.name[0]),
            end="",
        )
    num_atom += 1


#%%
def download_pdb(pdbid):
    """
    downloading pdb doesn't always succeed.
    """
    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    r = requests.get(url, allow_redirects=True)
    return r.content


def download_rcsb(filename):
    url = f"https://files.rcsb.org/download/{filename}"
    r = requests.get(url, allow_redirects=True)
    return r.content
