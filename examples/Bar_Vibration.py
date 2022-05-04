import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])
# Bar Vibration problem

# Lines for construction
C3=line(p0=(0.00,0.00,0), p1=(2.0*np.sqrt(2),0.00,0))
C4=line(p0=(0.00,2.0*np.sqrt(2),0), p1=(2.0*np.sqrt(2),2.0*np.sqrt(2),0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

S=sf.cent_stretch_bg(C3, C4, 0, np.array([10, 10, 0]), np.array([1.1,1.1,1.0]), np.array([1.0,1.0,1]), n=1)
#S=sf.nonlinear_parameterization(S) # Here the function nonlinear parameterization is called. This function preserves the order of the background as it is input!s

#print(S.control)
# Call foreground geometry generation with discritization and geometry parameters
G=ms.generate_unif_foreground(O, np.array([24.9998, 4.9998, 0.0]), np.array([50, 10, 1]))
G=sf.translate(G, 0.0001, 0.0001, 0)


# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G, S, 1)
ms.vis_background()
