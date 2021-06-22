import jnpy as jp
from jnpy import bio

pdb1 = bio.read_pdb('7joz.pdb')
pdb2 = bio.read_pdb('PW0441_WT D1R_dock.pdb')

aa1 = ['A',   'R',   'N',   'D',   'C',   'Q',   'E',   'G',   'H',   'I',   'L',   'K',   'M',   'F',   'P',   'S',   'T',   'W',   'Y',   'V',   'O',   'U',   'B',   'Z',   'X',   'J',   '*']
aa2 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'PYL', 'SEC', 'ASX', 'GLX', 'XAA', 'XLE', 'TERM']

def get_chain(model, chain):
  rs = []
  for c in model:
    if c.name == chain:
      for r in c:
        if r.name in aa2:
          rs.append(r)
  return rs

def seq_ind(seq):
  j = 0
  ind = []
  for c in seq:
    if c != '-':
      ind.append(j)
      j += 1
    else:
      ind.append(-1)
  return ind

def lig_rmsd(m1, m2):
  rs1 = get_chain(m1, 'R')
  rs2 = get_chain(m2, 'A')

  seq1 = ''.join(aa1[aa2.index(r.name)] for r in rs1)
  seq2 = ''.join(aa1[aa2.index(r.name)] for r in rs2)

  seq_aln1, seq_aln2, score = jp.align(seq1, seq2)

  common_indices = [i for i in range(len(seq_aln1)) if seq_aln1[i] != '-' and seq_aln2[i] != '-']
  ind1 = seq_ind(seq_aln1)
  ind2 = seq_ind(seq_aln2)

  rs1 = [rs1[ind1[i]] for i in common_indices]
  rs2 = [rs2[ind2[i]] for i in common_indices]

  sp, rmsd = jp.suppos(rs2, rs1, True)

  ligs1 = [r for c in m1 for r in c if r.name == 'VFP']
  ligs2 = [r for c in m2 for r in c if r.name == 'UNK']

  jp.sp_apply(sp, ligs2)
  rmsd = jp.rmsd(ligs2, ligs1, False, True)
  return rmsd

m1 = pdb1[0]
for m2 in pdb2:
  print(lig_rmsd(m1, m2))

