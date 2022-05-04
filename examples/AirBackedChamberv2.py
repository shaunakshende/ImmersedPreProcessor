import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])

padding = 0.5-0.305
# h = (0.305+padding)/165.0
h = (0.305+padding)/75.0
hf = (0.305+padding)/151.0
L = 0.305+padding
standoff = 0.076
# Lines for construction
C3=line(p0=(0.0,0.0,0), p1=(0.305+padding,0.0,0))
C4=line(p0=(0.0,0.305+padding,0), p1=(0.305+padding, 0.305+padding,0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

#S=sf.cent_stretch_bg(C3, C4, 0.305+1.0*padding, np.array([166, 166, 166]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)
S=sf.cent_stretch_bg(C3, C4, 0.305+1.0*padding, np.array([21, 21, 21]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)
# Add discontinuity at air/water interface.
#S.insert(2,5*h/(0.305+padding))

#print(S.control)
# Call foreground geometry generation with discritization and geometry parameters

# Plate (Homogenized Composite plate) - Material 0
#G0 = ms.generate_unif_PDforeground(O, np.array([0.305-0.00001, 0.305-0.00001, 0.013]), np.array([100, 100, 1]), 0)
#G0 = sf.translate(G0, 0.00001/2.0, 0.00001/2.0, 0.030+0.013/2.0)
G0 = ms.generate_unif_PDforeground(O, np.array([0.305, 0.305, 0.00126]), np.array([20, 20, 1]), 0)
G0 = sf.translate(G0, padding/2.0, padding/2.0, 0.2)

# Charge - Material 1 - RP-87 charge
#G1 = ms.generate_unif_foreground(O, np.array([0.004768, 0.004768, 0.002]), np.array([30, 30, 13]), 1)
#G1 = ms.Cylinder_Mask(G1, np.array([0.004768/2.0, 0.004768/2.0, 0.0]), 0.0024384)
#G1 = sf.translate(G1, 0.305/2-0.004768/2.0, 0.305/2-0.004768/2.0, 0.1+0.013+0.036) # 36 mm standoff distance

# Charge - Material 1 - RP-503 charge

G1 = ms.generate_unif_foreground(O, np.array([0.0051*2.0, 0.0051*2.0, 0.0047]), np.array([30, 30, 15]), 1)
G1 = sf.translate(G1, 0.305/2.0-0.0051+padding/2.0, 0.305/2.0-0.0051+padding/2.0,0.2+0.00126+standoff) # 76 mm standoff distance
G1 = ms.Cylinder_Mask(G1, np.array([L/2.0, L/2.0, 0]), 0.0049)

# Air - Material 2
#G2 = ms.generate_unif_foreground(O, np.array([0.305-0.00001, 0.305-0.00001, 0.030]), np.array([150, 150, 30]), 2)
#G2 = sf.translate(G2, 0.00001/2.0, 0.00001/2.0, 0.0)
# G2 = ms.generate_unif_foreground(O, np.array([0.305, 0.305, 10.0*h-h/4.0]), np.array([152, 152, 30]), 2)
# G2 = sf.translate(G2, padding/2.0, padding/2.0, 0.0)

#The Peridynamic solid has to be assembled last (input convention)
#G1 = sf.fg_superpose(G1, G2)
G1 = sf.fg_superpose(G1, G0)

# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G1, S, 4)
ms.save_PDGeometry(G0, 1, 'Plate', "NMS")
ms.vis_background()
