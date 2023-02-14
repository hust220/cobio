import plotly.graph_objects as go
import numpy as np
import math

import cobio as jp

p = jp.Pocket('1ec0.rec.pdb', '1ec0.lig.mol2', box=20, bin=5)

children = p.find_children()
ls1 = p.receptor_grid()
ls2 = [child.ligand_grid()[0] for child in children]

d = np.array(ls1+ls2)

fig = go.Figure(data=[go.Scatter3d(x=d.T[0], y=d.T[1], z=d.T[2], mode='markers')])
fig.show()


f68e8525a150c864992c4da4df0da7b3037933a55dafc0f3
