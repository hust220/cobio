import numpy as np
import math
import matplotlib.pyplot as plt
import jnpy as jp

p = jp.Pocket('1ec0.rec.pdb', '1ec0.lig.mol2', box=15, bin=1.5)

children = p.find_children()
print(len(children))
children2 = children[0].find_children()
ls1 = p.receptor_grid()
ls2 = [child.ligand_grid()[0] for child in children]

d1 = np.array(ls1)
d2 = np.array(ls2)

print(len(ls2))

