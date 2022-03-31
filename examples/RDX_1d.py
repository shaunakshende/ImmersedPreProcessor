import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])

# Lines for construction
C1=line(p0=(0.0,0.0,0.0), p1=(4,0.0,0.0))
C2=line(p0=(0.0,0.1,0.0), p1=(4,0.1,0.0))

C3=line(p0=(4,0.0,0.0), p1=(54.15,0.0,0.0))
C4=line(p0=(4,0.1,0.0), p1=(54.15,0.1,0.0))
#n=1 for background quadratic, 2 for cubic
Sc=sf.cent_stretch_bg(C1, C2, 0.0, np.array([4000+1, 1, 1]), np.array([0.25,1.1,1.0]), np.array([1.0,1.0,1.0]), n=1)
Sf=sf.cent_stretch_bg(C3, C4, 0.0, np.array([5000, 1, 1]), np.array([0.25,1.1,1.0]), np.array([1.0,1.0,1.0]), n=1)
S = sf.merge_bg(Sc, Sf, 0)

print(S.knots)
#Call foreground geometry generation with discritization and geometry parameters(0.0012-0.0012/3.0)
#0.0012/6.0
G=ms.generate_unif_foreground(O, np.array([0.005415, (0.1-0.1/3.0), 0.0]), np.array([500, 3, 1]))
G=sf.translate(G, 0.0, 0.1/6.0, 0)
ms.save_geometry(G, S, 1)
ms.vis_background()
#
