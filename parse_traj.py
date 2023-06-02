import sys
import math
import re
import os

class Atom:
  def __init__(self):
    self.name = ''
    self.num = 0
    self.coord = [0,0,0]
    self.is_het = False

  def __getitem__(self, i):
    return self.coord[i]

  def __setitem__(self, i, v):
    self.coord[i] = v

  def __iter__(self):
    return iter(self.coord)

  def distance(self, atom2):
    dx = self[0] - atom2[0]
    dy = self[1] - atom2[1]
    dz = self[2] - atom2[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)

class Residue:
  def __init__(self):
    self.name = ''
    self.num = 0
    self.atoms = []

  def __iter__(self):
    return iter(self.atoms)

  def __getitem__(self, key):
    if isinstance(key, int):
      return self.atoms[key]
    else:
      for atom in atoms:
        if atom.name == key:
          return atom
      print("Couldn't find atom '%s'" % key)
      write_residue(self)
      quit()

class Chain:
  def __init__(self):
    self.name = ''
    self.residues = []

  def __iter__(self):
    return iter(self.residues)

  def __getitem__(self, key):
    if isinstance(key, int):
      return super(Chain, self).__getitem__(key)
    else:
      for residue in self:
        if residue.num == key:
          return residue
      print("Couldn't find residue %s in chain %s" % (key, self.name))
      quit()

  @property
  def seq(self):
    return ''.join(res.name.strip()[0] for res in self.residues)

  @property
  def atoms(self):
    return [atom for residue in self.residues for atom in residue.atoms]

class Model:
  def __init__(self):
    self.num = 0
    self.chains = []

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
    return ''.join(res.name.strip()[0] for chain in self.chains for res in chain)

  @property
  def residues(self):
    return [residue for chain in self for residue in chain]

  @property
  def atoms(self):
    return [atom for chain in self for residue in chain for atom in residue]

class Pdb:
  def __init__(self):
    self.name = ''
    self.models = []

  def __iter__(self):
    return iter(self.models)

  def __getitem__(self, key):
    return self.models[key]

def read_pdb(filename):
  return PdbParser(filename).pdb

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
    self.is_het = line.startswith('HETATM')

def read_traj(filename):
  atom, residue, chain, model = Atom(), Residue(), Chain(), Model()

  def add_residue(l):
    if len(residue.atoms) > 0:
      residue.name = l.res_name
      residue.num = l.res_num
      chain.residues.append(residue)
      return Residue()
    else:
      return residue

  def add_chain(l):
    if len(chain.residues) > 0:
      chain.name = l.chain_name
      model.chains.append(chain)
      return Chain()
    else:
      return chain

  l, l_old = None, None
  for line in open(filename):
    if line.startswith('ATOM') or line.startswith('HETATM'):
      l = ParsedLine(line)

      if len(residue.atoms) > 0:
        # if residue changes
        if l.res_num != l_old.res_num or l.res_name != l_old.res_name or l.chain_name != l_old.chain_name:
          residue = add_residue(l_old)
          residue = Residue()
          # print(len(residue.atoms))
        # if chain changes
        if l.chain_name != l_old.chain_name:
          chain = add_chain(l_old)

      # if model changes
      atom.num = l.atom_num
      atom.name = l.atom_name
      atom.coord = [l.x, l.y, l.z]
      residue.atoms.append(atom)
      atom = Atom()

      l_old = l
      
    elif line.startswith('TER'):
      residue = add_residue(l_old)
      chain = add_chain(l_old)
    elif line.startswith('MODEL') or line.startswith('ENDMDL'):
      residue = add_residue(l_old)
      chain = add_chain(l_old)
      if len(model.chains) > 0:
        yield model
        model = Model()
    elif line.startswith('END'):
      residue = add_residue(l_old)
      chain = add_chain(l_old)
      if len(model.chains) > 0:
        return model
  residue = add_residue(l_old)
  chain = add_chain(l_old)
  if len(model.chains) > 0:
    return model
  else:
    return None

def mindist(r1, r2):
  mn = 9999
  for a1 in r1:
    for a2 in r2:
      d = a1.distance(a2)
      if d < mn:
        mn = d
  return mn

infile = sys.argv[1]
outfile = sys.argv[2]

l1 = [0]
l2 = list(range(36))
l_r1 = set()
l_af = set()

r1_atoms = ["P", "OP1", "OP2", "OP3", "HOP3", "O5'", \
"C5'", "H5'", "H5''", "C4'", "H4'", "C3'", "H3'", "O3'", \
"C2'", "H2'", "H2''", "C1'", "H1'", "O4'", \
"N1", "C2", "O2", "N3", "C4", "C5", "H5", "C6", "H6", "N4", "H41"]

traj = read_traj(infile)
distances = []
iframe = 0
f = open(outfile, 'w+')
for frame in traj:
  print('[*] Parsing model {0} ...'.format(iframe+1), flush=True)
#   calculate the distances based on the model
  dists = []
  rs = list(frame.residues)
  # print([len(r.atoms) for r in rs])
  if iframe == 0:
    names = [a.name for a in rs[0]]
    # print(names)
    l_r1 = set(names.index(name) for name in r1_atoms)
    l_af = set(range(len(names))) - l_r1
  rs1 = [[rs[0][i] for i in l_r1]]
  rs2 = [[rs[0][i] for i in l_af]] + [rs[i] for i in range(1,36)]
  # ir = 0
  # for r in rs2:
  #   for a in r:
  #     print(ir+1, a.name, a.num, a.coord)
  #   ir += 1
  for r1 in rs1:
    for r2 in rs2:
      d = mindist(r1, r2)
      dists.append(d)
  # print(dists)
  f.write('{0}\n'.format(' '.join(str(d) for d in dists)))

  distances.append(dists)
  iframe += 1

f.close()
#     save_data()
#     plot_figure()