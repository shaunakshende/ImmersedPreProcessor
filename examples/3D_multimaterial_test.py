import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])


# Lines for construction
C3=line(p0=(0.0,0.00,0), p1=(0.305,0.00,0))
C4=line(p0=(0.0,0.305,0), p1=(0.305, 0.305,0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

S=sf.cent_stretch_bg(C3, C4, 0.610, np.array([101, 101, 201]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)

#print(S.control)
# Call foreground geometry generation with discritization and geometry parameters

# Plate (Homogenized Composite plate) - Material 0
G0 = ms.generate_unif_PDforeground(O, np.array([0.305-0.00001, 0.305-0.00001, 0.013]), np.array([10, 10, 1]), 0)
G0 = sf.translate(G0, 0.00001/2.0, 0.00001/2.0, 0.030+0.013/2.0)

# Charge - Material 1 - RP-87 charge
#G1 = ms.generate_unif_foreground(O, np.array([0.004768, 0.004768, 0.002]), np.array([30, 30, 13]), 1)
#G1 = ms.Cylinder_Mask(G1, np.array([0.004768/2.0, 0.004768/2.0, 0.0]), 0.0024384)
#G1 = sf.translate(G1, 0.305/2-0.004768/2.0, 0.305/2-0.004768/2.0, 0.1+0.013+0.036) # 36 mm standoff distance 

# Charge - Material 1 - RP-503 charge
G1 = ms.generate_unif_foreground(O, np.array([0.0051*2.0, 0.0051*2.0, 0.0047]), np.array([5, 5, 5]), 1)
G1 = ms.Cylinder_Mask(G1, np.array([0.0051, 0.0051, 0]), 0.0049)
G1 = sf.translate(G1, 0.305/2.0-0.0049/2.0, 0.305/2.0-0.0049/2.0, 0.030+0.013+0.01) # 10 mm standoff distance

# Air - Material 2
G2 = ms.generate_unif_foreground(O, np.array([0.305-0.00001, 0.305-0.00001, 0.030]), np.array([15, 15, 3]), 2)
G2 = sf.translate(G2, 0.00001/2.0, 0.00001/2.0, 0.0)

#The Peridynamic solid has to be assembled last (input convention)
G1 = sf.fg_superpose(G1, G2)
G1 = sf.fg_superpose(G1, G0)

# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G1, S, 192)
ms.save_PDGeometry(G0, 1, 'Plate')
ms.vis_background()
