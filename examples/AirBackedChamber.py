import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])


# Lines for construction
C3=line(p0=(0.0,0.00,0), p1=(0.305,0.00,0))
C4=line(p0=(0.0,0.6,0), p1=(0.305, 0.6,0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

S=sf.cent_stretch_bg(C3, C4, 0.305, np.array([11, 16, 11]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)
#S=sf.nonlinear_parameterization(S) # Here the function nonlinear parameterization is called. This function preserves the order of the background as it is input!s

#print(S.control)
# Call foreground geometry generation with discritization and geometry parameters

# Plate (Steel) - Material 0
# Because of the way Peridigm assigns points to processors, the total number of particles must be divisible
# by the number of processors!!!!! I.e. the number of particles on each processor are the SAME or
# the PD_ID used to map quantites from PD to IGA-Blast will be incorrect!
G0 = ms.generate_unif_PDforeground(O, np.array([0.305/2.0, 0.305, 0.013]), np.array([10, 30, 1]), 0)
G0 = sf.translate(G0, 0.0, 0.5*(0.6-0.305), 0.1)

# Charge - Material 1
G1 = ms.generate_unif_foreground(O, np.array([0.004768, 0.004768, 0.002]), np.array([10, 10, 10]), 1)
G1 = ms.Cylinder_Mask(G1, np.array([0.0, 0.004768/2.0, 0.0]), 0.0024384)
G1 = sf.translate(G1, 0.0, 0.6/2.0-0.004768/2.0, 0.1+0.013+0.036) # 100 mm standoff distance

# Air - Material 2
#G2 = ms.generate_unif_foreground(O, np.array([0.305/2.0, 0.305, 0.1]), np.array([10, 20, 7]), 2)
#G2 = sf.translate(G2, 0.0, 0.5*(0.6-0.305), 0.0)


G1 = sf.fg_superpose(G1, G2)
#The Peridynamic solid has to be assembled last (input convention)
G1 = sf.fg_superpose(G1, G0)

# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G1, S, 3)
ms.save_PDGeometry(G0, 1, 'Plate', "NMS")
ms.vis_background()
