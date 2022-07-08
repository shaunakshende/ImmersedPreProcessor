import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])

h = 2.4/119.0;
padding_x = (2.4-1.2)/2.0
padding_z = 120.0/2.0*h-0.32

# Lines for construction
C3=line(p0=(0.0,0.0,0), p1=(2.4,0.0,0))
C4=line(p0=(0.0,2.4,0), p1=(2.4, 2.4, 0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

S=sf.cent_stretch_bg(C3, C4, 2.4, np.array([120, 120, 120]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)

# Material 0 = Concrete
# G0 = ms.generate_unif_PDforeground(O, np.array([1.2, 1.2, 0.32]), np.array([150, 150, 40]), 0)
# G0 = sf.translate(G0, padding_x, padding_x, padding_z)
G0 = ms.generate_unif_PDforeground(O, np.array([0.6, 0.6, 0.32]), np.array([120/2, 120/2, 34]), 0)

#G0 = ms.generate_unif_PDforeground(O, np.array([1.2, 1.2, 0.32]), np.array([12, 12, 3]), 0)
G0 = sf.translate(G0, padding_x*0.0, padding_x*0.0, padding_z)

# # Charge - Material 1 - TNT charge
G1 = ms.generate_FGCone_cylindrical2(O, 0.103/2.0, 0.0085, 0.001, 0.075, [40, 200/4, 50], 1, 2, 1, angle = np.pi/2.0)
G1 = sf.translate(G1, 0.0, 0.0, 0.32+padding_z+h/4.0)

# # Material 2 : CompB (mix of RDX and TNT)
G2 = ms.generate_FGCone_cylindrical2(O, 0.113/2.0, 0.009325, 0.001, 0.08228, [44, 220/4, 55], 1, 2, 2, angle = np.pi/2.0)
G2 = ms.Conic_Mask2(G2, O, 0.075, 0.103/2.0, 0.0085, 2, 0)
G2 = sf.translate(G2, 0.0, 0.0, 0.32+padding_z+h/4.0)
#
G1 = sf.fg_superpose(G1, G2)
# #The Peridynamic solid has to be assembled last (input convention)
G1 = sf.fg_superpose(G1, G0)
# # Save geometry for visualization, and background as .dat files.

ms.save_geometry(G1, S, 800)
ms.save_PDGeometry(G0, 1, 'Block', "mmNS")
ms.vis_background()
