import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])

h = 0.3/15.0;
padding_x = (1.8-1.2)/2.0
padding_z = 90.0/2.0*h-0.32

# Lines for construction
C3=line(p0=(0.0,0.0,0), p1=(0.3,0.0,0))
C4=line(p0=(0.0,0.3,0), p1=(0.3, 0.3, 0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

S=sf.cent_stretch_bg(C3, C4, 0.3, np.array([16, 16, 16]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)

G0 = ms.generate_unif_PDforeground(O, np.array([0.305, 0.305, 0.00126]), np.array([10, 100, 20]), 0)
#G0 = sf.translate(G0, padding/2.0, padding/2.0, 0.00001*h)
# Charge - Material 1 - RDX/ls
# TNT/CompB shape charge
#[layers, circumferential points, layers height]
G1 = ms.generate_FGCone_cylindrical2(O, 0.103/2.0, 0.0085, 0.001, 0.075, [40, 200, 50], 1, 2, 1)
G1 = sf.translate(G1, 0.9, 0.9, 0.0)#0.32+padding_z+h/4.0)

#G2 = ms.generate_FGCone_cylindrical2(O, 0.113/2.0, 0.009325, 0.001, 0.08228, [44, 220, 55], 1, 2, 2)
#G2 = ms.Conic_Mask2(G2, O, 0.075, 0.103/2.0, 0.0085, 2, 0)
#G2 = sf.translate(G2, 0.9, 0.9, 0.0)#0.32+padding_z+h/4.0)
#The Peridynamic solid has to be assembled last (input convention)
#G1 = sf.fg_superpose(G1, G2)
# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G1, S, 5*80)
ms.save_PDGeometry(G0, 1, 'Plate', "NMS")
ms.vis_background()
