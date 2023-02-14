import cobio as jp

f = open('bb.pdb', 'r')
pdb = jp.Pdb(f)
f.close()
print(pdb.seq())
#print(pdb.contacts(cutoff=5))
