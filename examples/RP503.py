import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])


# Lines for construction
C3=line(p0=(0.0,0.00,0), p1=(0.31,0.00,0))
C4=line(p0=(0.0,0.31,0), p1=(0.31, 0.31,0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

#S=sf.cent_stretch_bg(C3, C4, 1.0, np.array([100, 100, 100]), np.array([0,0,0]), np.array([1.16,1.16,1.55]), n=1)

S=sf.cent_stretch_bg(C3, C4, 0.31, np.array([191, 191, 191]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)

# Call foreground geometry generation with discritization and geometry parameters
G = ms.generate_unif_foreground(O, np.array([0.0051, 0.0051, 0.0047*0.5]), np.array([70, 70, 70]))

G = ms.Cylinder_Mask(G, np.array([0,0,0]), 0.0049)

G = sf.translate(G, 0.00000, 0.00000, 0.00000)
# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G, S, 4*96)
ms.vis_background()
