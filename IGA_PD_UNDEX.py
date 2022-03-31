import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])


# Lines for construction
C3=line(p0=(0.0,0.00,0), p1=(0.4,0.00,0))
C4=line(p0=(0.0,0.4,0), p1=(0.4, 0.4,0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

S=sf.cent_stretch_bg(C3, C4, 0.0, np.array([80, 80, 1]), np.array([0.0,0.0,0]), np.array([1.2,1.2,1.0]), n=1)
#S=sf.nonlinear_parameterization(S) # Here the function nonlinear parameterization is called. This function preserves the order of the background as it is input!s

#print(S.control)
# Call foreground geometry generation with discritization and geometry parameters
G = ms.generate_unif_foreground(O, np.array([0.2, 0.1, 0.0]), np.array([15, 20, 1]))
#G = sf.cent_stretch_fg(np.array([0.04, 0.05, 0.0]),np.array([1.2,1.2,1.0]), G)
# Spherical mask on the unform rectallinear domain
#G = ms.Sphere_Mask(G, [0.0, 0.0, 0.0], 0.03)

G = sf.translate(G, 0.1, 0.15, 0.00000)

# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G, S, 1)
ms.vis_background()
