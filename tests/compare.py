import cobio as jp
from cobio import bio

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

  as1 = [a for r in ligs1 for a in r]
  as2 = [a for r in ligs2 for a in r]

  mp = jp.map_atoms(as1, as2)
  f1 = open('lig1.pdb', 'w+')
  f2 = open('lig2.pdb', 'w+')
  for i in range(len(as1)):
    a1 = as1[i]
    a2 = as2[mp[i]]
    f1.write("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n" % \
      (i+1, a1.name[0]+str(i+1), 'VFP', 'A', 1, a1[0], a1[1], a1[2] , 1.00 , 0.00, a1.name[0]))
    f2.write("ATOM%7i  %-4s%3s%2s%4i%12.3lf%8.3lf%8.3lf%6.2f%6.2f%12c  \n" % \
      (i+1, a2.name[0]+str(i+1), 'UNK', 'A', 1, a2[0], a2[1], a2[2] , 1.00 , 0.00, a2.name[0]))
  f1.close()
  f2.close()

  rmsd = jp.rmsd(ligs2, ligs1, False, True)
  return rmsd

m1 = pdb1[0]
for m2 in pdb2:
  print(lig_rmsd(m1, m2))

