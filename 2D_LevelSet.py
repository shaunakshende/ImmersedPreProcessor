import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])

# Lines for construction
C3=line(p0=(0.0,0.0,0), p1=(0.5,0.0,0))
C4=line(p0=(0.0,0.5,0), p1=(0.5, 0.5, 0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

S=sf.cent_stretch_bg(C3, C4, 0, np.array([150, 150, 1]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)

# Material 0 = currently dummy material
# G0 = ms.generate_unif_PDforeground(O, np.array([1.2, 1.2, 0.32]), np.array([150, 150, 40]), 0)
# G0 = sf.translate(G0, padding_x, padding_x, padding_z)
G0 = ms.generate_unif_PDforeground(O, np.array([1.2, 1.2, 0.32]), np.array([10, 10, 1]), 0)

#The Peridynamic solid has to be assembled last (input convention)
ms.save_geometry(G0, S, 1)
ms.vis_background()
