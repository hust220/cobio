import jnpy as jp
p = jp.Pocket('1ec0.rec.pdb', '1ec0.lig.mol2', box=10, bin=0.25)
print(p.ligand_grid())
print(p.receptor_grid())
children = p.find_children()

print(children[0].ligand_grid())

children[0].find_children()

#print(children[0].ligand_grid())
