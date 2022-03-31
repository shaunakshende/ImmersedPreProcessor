import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])


# Lines for construction
C3=line(p0=(0.0,0.00,0), p1=(2.0,0.00,0))
C4=line(p0=(0.0,2.0,0), p1=(2.0, 2.0,0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

S=sf.cent_stretch_bg(C3, C4, 2.0, np.array([25, 25, 25]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)
#S=sf.nonlinear_parameterization(S) # Here the function nonlinear parameterization is called. This function preserves the order of the background as it is input!s

#print(S.control)
# Call foreground geometry generation with discritization and geometry parameters
G = ms.generate_unif_foreground(O, np.array([0.161, 0.161, 0.161]), np.array([15, 15, 15]))

#G2 = ms.generate_unif_foreground(O, np.array([0.161, 0.161, 0.161]), np.array([10, 10, 10]))
#Apply some foreground stretching to account for the background stretching as well.

#G = sf.cent_stretch_fg(O,np.array([1.,1.35,1.35]), G)

# Spherical mask on the unform cubic domain
G = ms.Sphere_Mask(G, O, 0.16)
#G2 = ms.Sphere_Mask(G, O, 0.16)
#G2 = sf.translate(G2, 1.5, 1.5, 1.5)

#G3 = sf.fg_superpose(G, G2)

# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G, S, 3)
ms.vis_background()
