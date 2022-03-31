import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np
# Specify origin
O=np.array([0,0,0])

# Lines for construction
C3=line(p0=(0.00,0.00,0), p1=(0.003,0.00,0))
C4=line(p0=(0.00,0.01,0), p1=(0.003,0.01,0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
S=sf.cent_stretch_bg(C3, C4, 0, np.array([98,198,1]), np.array([0.003,0.005,0]), np.array([1.0,1.0,1]), n=2)

# Call foreground geometry generation with discritization and geometry parameters
G=ms.generate_unif_foreground(O, np.array([0.008,0.008,0.0]), np.array([50,50,1]))
G=ms.subt_circular_domain(G, np.array([0.004,0.004,0]), 0.0005)
ms.vis_foreground2d(G)
#Biaxial stretching
G=sf.lin_stretch(0, 0.002, 1.0, G)
G=sf.lin_stretch(1, 0.004, 1.0, G)
G=sf.translate(G,0.00099999,0.001,0)

# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G, S, 1)
ms.vis_background()
