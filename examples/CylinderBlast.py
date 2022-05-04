import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])


G = ms.generate_unif_foreground(O, np.array([0.05, 0.05, 0.001]), np.array([30, 30, 1]))
G = sf.translate(G, 0.00000, 0.00000, 0.00000)

ms.save_PDGeometry(G, 1, 'CylinderBlast')

