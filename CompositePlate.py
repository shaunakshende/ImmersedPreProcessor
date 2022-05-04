import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np

# Specify origin
O=np.array([0,0,0])

padding = 0.5-0.305
h = (0.305+padding)/75.0
hf = (0.305+padding)/151.0
L = 0.305+padding
standoff = 0.152
# Lines for construction
C3=line(p0=(0.0,0.0,0), p1=(0.305+padding,0.0,0))
C4=line(p0=(0.0,0.305+padding,0), p1=(0.305+padding, 0.305+padding,0))

# Call geometry generation functions for the background with optional argument n which indicates elevation of background elements
# n=1 for background quadratic, 2 for cubic

#S=sf.cent_stretch_bg(C3, C4, 0.305+1.0*padding, np.array([101, 101, 101]), np.array([0,0,0]), np.array([1.0,1.0,1.0]), n=1)
S = ms.generate_unif_background(C3, C4, 0.305+1.0*padding, np.array([201, 201, 201]));
# Call foreground geometry generation with discritization and geometry parameters
# Plate (Homogenized Composite plate) - Material 0
G0 = ms.generate_unif_PDforeground(O, np.array([0.305, 0.305, 0.00126]), np.array([60, 60, 1]), 0)
#G0 = ms.generate_unif_PDforeground(O, np.array([0.305, 0.305, 0.00126]), np.array([366, 366, 1]), 0)
G0 = sf.translate(G0, padding/2.0, padding/2.0, 0)
print(np.sum(G0.vols)-0.305*0.305*0.00126)

# Charge - Material 1 - RP-503 charge
G1 = ms.generate_unif_foreground(O, np.array([0.0051*2.0, 0.0051*2.0, 0.004651320879015]), np.array([71, 71, 33]), 1)
G1 = sf.translate(G1, 0.305/2.0-0.0051+padding/2.0, 0.305/2.0-0.0051+padding/2.0, 0.150-(0.004651320879015)/2.0) # 152 mm standoff distance
G1 = ms.Cylinder_Mask(G1, np.array([L/2.0, L/2.0, 0]), 0.0049)
vol = G1.vols;

print(abs(np.sum(vol)-np.pi*0.0049**2*0.004651320879015))
print(abs(np.sum(vol)*1770.0*1000000))
print(np.pi*0.0049**2*0.004651320879015*1770.0*1000000)

# Air - Material 2
#G2 = ms.generate_unif_foreground(O, np.array([0.305-0.00001, 0.305-0.00001, 0.030]), np.array([150, 150, 30]), 2)
# #G2 = sf.translate(G2, 0.00001/2.0, 0.00001/2.0, 0.0)
# G2 = ms.generate_unif_foreground(O, np.array([0.305, 0.305, 10.0*h-h/4.0]), np.array([152, 152, 15]), 2)
# G2 = sf.translate(G2, padding/2.0, padding/2.0, 0.0)

#The Peridynamic solid has to be assembled last (input convention)
#G1 = sf.fg_superpose(G1, G2)
G1 = sf.fg_superpose(G1, G0)

# Save geometry for visualization, and background as .dat files.
ms.save_geometry(G1, S, 240)

# save PD geometry : Foreground object, # processors, name, units (NMS = newton meter second)
ms.save_PDGeometry(G0, 1, 'Plate', "NMS")
ms.vis_background()
