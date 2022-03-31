from igakit.cad import *
import numpy as np
from igakit.io import PetIGA,VTK
from numpy import linspace
Geo=PetIGA().read("./Geometry.dat")
print(Geo.degree)
