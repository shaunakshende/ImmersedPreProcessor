import mesh_functions as ms
import stretching_functions as sf
from igakit.cad import *
import numpy as np
O=np.array([0,0,0])


C3=line(p0=(0.0000,0.00,0), p1=(0.003,0.00,0))
C4=line(p0=(0.0000,0.01,0), p1=(0.003,0.01,0))
S=ms.generate_unif_background(C3, C4, 0, np.array([52,167,1]))



G=ms.generate_unif_foreground(O, np.array([0.002,0.008,0.0]), np.array([100,400,1]))
G=lin_stretch(1, 0, 1.01, 0)

ms.save_geometry(G.coor, G.vols, S, 1)


ms.vis_background()

